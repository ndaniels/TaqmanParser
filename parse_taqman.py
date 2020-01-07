#!/usr/bin/env python3

import sys
import os
import re

"""
Requires 1 command line argument: input filename.
Takes optional second command line argument: output filename. Missing this, it writes to stdout.
This program parses the Taqman output file in order to do SNP calls per Rachel F. Daniels' specifications.
It requires the following fields be present:
Under the [Sample Setup] section, the fields:
    "Sample Name", “SNP Assay Name”,  "Allele1 Name”, and “Allele2 Name”
must be present.
Under the [Results] section, the fields:
    “Sample Name”, “SNP Assay Name”, "Allele1 Crt”, and "Allele2 Crt”
must be present.
The combination of Sample Name and SNP Assay Name is a unique identifier and used as a
foreign key to join the two data sections (Sample Setup x Results).

The parsing works as a state machine, where states are in_setup and in_results.
The number of (tab-delimited) fields in each section is determined by the (tab-delimited) header.
If the number of fields in a line differs, a runtime error is raised.
If a blank line or EOF is encountered while in in_results or in_setup, that state is exited.
When in_results and in_setup are both False, but got_setup_header and got_results_header are both True,
then data have all been read and can be output.

[Sample Setup] —> "Sample Name", “SNP Assay Name”,  "Allele1 Name”, and “Allele2 Name” will give you the list of barcode positions (A1, B1, A2, B2… A12, B12) and their two possibly allele values (ATCG or ACTG)

[Results] —> “Sample Name”, “SNP Assay Name”, "Allele1 Crt”, and "Allele2 Crt”. SNP Assay Name is the key between them, as the name is identical between them. 
"""


"""
Given a taqman file name,
Parse the taqman file, and return two dictionaries.
For both dictionaries, the key is sample name concatenated with assay name.
Values in the setup dictionary are tuples of strings: (sample name, assay name, allele 1 name, allele 2 name)
Values in the results dictionary are tuples of strings: (allele 1 Ctr value, allele 2 Ctr value) 
"""
def parse_file(fname):
    setup = {}
    results = {}
    fname = sys.argv[1]
    assert(os.path.isfile(fname))
    with open(fname) as fh: 
        in_setup = False
        got_setup_header = False
        setup_header = []
        in_results = False
        got_results_header = False
        results_header = []
        for lnum, line in enumerate(fh,1):
            """
            Without loss of generality, handle Sample Setup first.
            """
            if re.match(r'.*\[Sample Setup\]', line):
                in_setup = True
                continue
            if in_setup and not got_setup_header:
                setup_header = re.split(r'\t', line)
                if setup_header:
                    got_setup_header = True
                    continue
            # Check if we should exit in_setup state
            if in_setup and len(line) <= 1: in_setup = False
            if in_setup and got_setup_header:
                fields = re.split(r'\t', line)
                assert (len(fields) == len(setup_header)), \
                    f"In Sample Setup, header has {len(setup_header)} fields but row has {len(fields)} at line {lnum}:\n{line}"
                sname, aname, a1name, a2name = '','','',''
                for h,f in zip(setup_header, fields):
                    if h == "Sample Name": sname = f
                    if h == "SNP Assay Name": aname = f
                    if h == "Allele1 Name": a1name = f
                    if h == "Allele2 Name": a2name = f
                assert (len(sname) >= 0 and len(aname) >= 0 and len(a1name) >= 0 and len(a2name) >= 0), f"Empty field on line {lnum}"
                assert sname+aname not in setup, f"Setup has multiple entries for {sname} {aname}"
                setup[sname+aname] = (sname, aname, a1name, a2name)
                
            """
            Results section
            """
            if re.match(r'.*\[Results\]', line):
                in_results = True
                continue
            if in_results and not got_results_header:
                # results section
                results_header = re.split(r'\t', line)
                if results_header:
                    got_results_header = True
                    continue
            # Check if we should exit in_results state
            if in_results and len(line) <= 1: in_results = False
            if in_results and got_results_header:
                fields = re.split(r'\t', line)
                assert len(fields) == len(results_header), \
                    f"In Results, header has {len(results_header)} fields but row has {len(fields)} at line {lnum}"
                # “Sample Name”, “SNP Assay Name”, "Allele1 Crt”, and "Allele2 Crt”
                sname, aname, a1crt, a2crt = '','','',''
                for h,f in zip(results_header, fields):
                    if h == "Sample Name": sname = f
                    if h == "SNP Assay Name": aname = f
                    if h == "Allele1 Crt": a1crt = f
                    if h == "Allele2 Crt": a2crt = f
                assert (len(sname) >= 0 and len(aname) >= 0 and len(a1crt) >= 0 and len(a2crt) >= 0), f"Empty field on line {lnum}"
                assert sname+aname not in results, f"Results has multiple entries for {sname} {aname}"
                results[sname+aname] = (a1crt, a2crt)
        """
        At this point we have read every line of the file
        """
        in_results = False
        in_setup = False
        assert len(results) == len(setup), \
            f"Setup and Results have differing numbers of rows! Setup: {len(setup)}; Results: {len(results)}"
        assert got_setup_header, "Did not find Setup section!"
        assert got_results_header, "Did not find Results section!"
    return (setup, results)

"""
Handle SNP-calling logic
If both Crt values are Undetermined, it's missing data (X).
If one is Undetermined, the other one is the call.
If neither is Undetermined, then if they are closer than 6 units apart, it's heterozygous (N).
Otherwise, the smaller is the call.
Returns a dictionary of calls, where the key is sample name concatenated with assay name.
"""
def call_snps(setup, results):
    threshold = 6.0
    het = 'N'
    missing = 'X' 
    undet = "Undetermined"
    calls = {}
    for key in setup.keys():
        sname, aname, a1name, a2name = setup[key]
        assert key in results, f"Results missing entry for {sname} {aname}"
        a1crt, a2crt = results[key]
        if a1crt == undet and a2crt == undet:
            calls[key] = missing
        elif a1crt == undet:
            calls[key] = a2name
        elif a2crt == undet:
            calls[key] = a1name
        else:
            a1f, a2f = float(a1crt), float(a2crt)
            if abs(a1f-a2f) < threshold:
                calls[key] = het
            elif a1f < a2f:
                calls[key] = a1name
            else:
                calls[key] = a2name
    return calls

"""
Format the results so that every sample is associated with its barcode (the calls in assay_sort order)
Returns a list of (samplename, barcode) tuples, sorted by samplename.
"""
def format_results(setup, calls):
    assays = assay_sort(setup)
    samples = get_samples(setup)
    sample_calls = []
    for sample in samples:
        barcode = ""
        for assay in assays:
            barcode += calls[sample+assay]
        sample_calls.append((sample, barcode))
    sample_calls.sort(key = lambda t: t[0])
    return sample_calls



"""
Given the setup dictionary, return the list of unique samples seen.
"""
def get_samples(setup):
    samples = set()
    for (k,v) in setup.items():
        sample, _, _, _ = v
        samples.add(sample)
    return list(samples)


"""
Look at the set of all assays seen in Sample Setup, and sort them in a manner
that is specific to the 384-well plates used for the P. falciparum barcode.
(A1, B1, A2, B2, A3, B3, ... A12, B12)
"""
def assay_sort(setup):
    seen_assays = set()
    for k,v in setup.items():
        _, assay, _, _ = v
        seen_assays.add(assay)
    """
    Sort the assays so they are sorted first by the digits and then by the character in the first position
    """
    assays = [a[1:]+a[0] for a in list(seen_assays)]
    assays.sort(key = lambda s: [int(t) if t.isdigit() else t.lower() for t in re.split('(\d+)', s)])
    """
    And restore the proper assay names
    """
    return [a[-1]+a[:-1] for a in assays]

def main():
    assert len(sys.argv) > 1, f"Usage: {sys.argv[0]} <input> [output]"
    infilename = sys.argv[1]
    outfilename = None
    if len(sys.argv) == 3:
        outfilename = sys.argv[2]
    setup, results = parse_file(outfilename)
    calls = call_snps(setup, results)
    formatted_results = format_results(setup, calls)
    output = ""
    for sample, barcode in formatted_results:
        output += f"{sample}\t{barcode}\n"
    if outfilename:
        with open(outfilename, 'w') as h:
            h.write(output)
    else:
        sys.stdout.write(output)

main()
