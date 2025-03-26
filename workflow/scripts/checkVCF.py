#!/usr/bin/env python3
"""
checkVCF.py -- Check validity of VCF files

Original Author: Xiaowei Zhan (zhanxw@umich.edu), Dajiang Liu (dajiang@umich.edu)
Original Source: https://github.com/zhanxw/checkVCF (version 1.4, 20140115)
License: No Licence specified, but repo was public.
Adapted by: Dazcam (20250321)
Further Modified: (20250324) to output rsIDs to --exclude file instead of chrom:pos

This Python 3 version is adapted from the original Python 2 script available at the above GitHub repository.
The original script validates VCF files by checking for duplicates, non-SNPs, reference mismatches,
invalid genotypes, monomorphic sites, and allele frequencies > 0.5, reporting issues to separate files.

Modifications in this version:
1. **Python 3 Compatibility**:
   - Updated syntax (e.g., print statements to functions with `file=`).
   - Handled bytes/string distinctions (e.g., `myopen` uses text mode for VCF/.fai).
   - Removed redundant Python 2 compatibility code (e.g., custom `all`, `os.SEEK_SET`).

2. **Command-Line Options**:
   - Added `-v` to explicitly specify the VCF file (original used positional argument).
   - Added `--exclude` to write all failing SNPs (originally chrom:pos, now rsIDs) to a single file.

3. **Exclusion Reporting**:
   - SNPs failing any check (non-SNPs, duplicates, reference mismatches, invalid genotypes,
     monomorphic sites, AF > 0.5) are written to the `--exclude` file as rsIDs if provided.
   - Original separate output files (e.g., `.check.dup`, `.check.ref`) are retained.

4. **Behavior**:
   - The script does NOT modify the input VCF; it only reports issues.
   - Added explicit UTF-8 decoding for consistency in byte/string handling.
   - Modified to output rsIDs instead of chrom:pos to `--exclude` file (20250324).

Contact for issues with original script: zhanxw@umich.edu or dajiang@umich.edu.
"""

import sys
import os
import logging
import getopt

VERSION = "version 2.0 (20250324)"  # Updated date to reflect this change

# Convenient functions
def myopen(fn):
    import gzip
    f = gzip.open(fn.encode(), 'rb')
    try:
        f.read(2)
        f.close()
        return gzip.open(fn.encode(), 'rt')  # Text mode for VCF/.fai
    except:
        f.close()
        return open(fn, 'r')  # Text mode

def checkGTformat(gt):
    if gt == '.':
        return True
    genos = gt.replace('|', '/').split('/')
    if len(genos) != 1 and len(genos) != 2:
        return False
    for g in genos:
        if g == '.':
            continue
        if g.isdigit():
            continue
        return False
    return True

def getGeno(gt):
    if gt.find('.') >= 0:
        return -1
    genos = gt.replace('|', '/').split('/')
    return sum([int(g) for g in genos])

# GenomeSequence class for .fai file handling
class GenomeSequence:
    def __init__(self):
        self.fn = ''
        self.data = {}
        self.fileHandle = None

    def open(self, fn):
        if not os.path.exists(fn + '.fai'):
            print(f"Cannot find .fai index file for {fn}, consider creating it with 'samtools faidx {fn}'", file=sys.stderr)
            return False
        self.fn = fn
        try:
            self.fileHandle = open(self.fn, 'rb')
            for ln in myopen(fn + '.fai'):
                fd = ln.strip().split()
                if fd[0][:3].upper() == "CHR":
                    fd[0] = fd[0][3:]
                key = fd[0]
                val = list(map(int, fd[1:]))
                self.data[key] = val
                self.data['CHR' + key] = val
            return True
        except Exception as e:
            print(e, file=sys.stderr)
            return False

    def getBase(self, chrom, pos):
        chrom = chrom.upper()
        if chrom not in self.data:
            return None
        size, loc, basesPerLine, bytePerLine = self.data[chrom]
        if pos < 0 or pos >= size:
            return None
        lineNo = pos // basesPerLine
        remainder = pos % basesPerLine
        filePos = loc + lineNo * bytePerLine + remainder
        self.fileHandle.seek(filePos, os.SEEK_SET)
        return self.fileHandle.read(1)

    def close(self):
        if self.fileHandle:
            self.fileHandle.close()

# Logging and UI functions
def actionItem(logger=sys.stderr):
    print("---------------     ACTION ITEM     ---------------", file=logger)

def usage():
    print("Usage: ")
    print(f"{sys.argv[0]} -r ref.fa -o prefix -v input.vcf [--exclude exclude.list]: check VCF for strand")

class Logger:
    def __init__(self, fn):
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s %(message)s',
                            datefmt='%m-%d %H:%M',
                            filename=fn,
                            filemode='w')
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        formatter = logging.Formatter('%(message)s')
        console.setFormatter(formatter)
        logging.getLogger('').addHandler(console)

    def write(self, msg):
        if msg == '\n':
            return
        logging.info(msg)

def banner(logger=sys.stderr):
    print("checkVCF.py -- check validity of VCF file for meta-analysis", file=logger)
    print(VERSION, file=logger)
    print("contact zhanxw@umich.edu or dajiang@umich.edu for problems.", file=logger)

# Main execution
if __name__ == '__main__':
    try:
        optlist, args = getopt.getopt(sys.argv[1:], 'r:o:v:', ['exclude='])
        optlist = dict(optlist)
        refFile = optlist['-r']
        outPrefix = optlist['-o']
        vcfFile = optlist.get('-v', args[0] if args else None)
        excludeFile = optlist.get('--exclude')
        if not vcfFile:
            raise ValueError("VCF file must be provided with -v or as a positional argument")
    except Exception as e:
        usage()
        print(e, file=sys.stderr)
        sys.exit(1)

    gs = GenomeSequence()
    if not gs.open(refFile):
        print("Cannot open reference genome, exiting...", file=sys.stderr)
        sys.exit(1)

    logger = Logger(outPrefix + '.check.log')
    banner(logger)

    # Open output files
    fDup = open(outPrefix + '.check.dup', 'w')
    fRef = open(outPrefix + '.check.ref', 'w')
    fNonSnp = open(outPrefix + '.check.nonSnp', 'w')
    fMono = open(outPrefix + '.check.mono', 'w')
    fGeno = open(outPrefix + '.check.geno', 'w')
    fAF = open(outPrefix + '.check.af', 'w')
    fExclude = open(excludeFile, 'w') if excludeFile else None

    nRef, nGeno, nAF = 0, 0, 0
    nNonSnp = 0
    snpSite = set()
    nDupSite = 0
    nMono = 0

    nField = -1
    ACGT = set(['A', 'C', 'G', 'T'])
    ACGTM = set(['A', 'C', 'G', 'T', '.'])
    lineNo = -1

    vcfHeaderLine = 0
    vcfSiteLine = 0
    vcfSample = 0
    chrWarningGiven = False

    print(f"Python version is [ {'.'.join(map(str, sys.version_info)).strip()} ]", file=logger)
    print(f"Begin checking vcfFile [ {vcfFile} ]", file=logger)

    try:
        prevChrom, prevPos = None, None
        for lineNo, ln in enumerate(myopen(vcfFile)):
            if lineNo % 10000 == 0 and lineNo != 0:
                print(f"[ {lineNo} ] lines processed \r", file=logger, end="")
            if not ln or ln.startswith('##'):
                vcfHeaderLine += 1
                continue
            if ln.startswith('#CHROM'):
                vcfHeaderLine += 1
                fd = ln.strip().split()
                nField = len(fd)
                if len(set(fd[9:])) != len(fd[9:]):
                    actionItem(logger)
                    print("Your VCF file has duplicated sample IDs, please fix them and re-run checkVCF.py", file=logger)
                    sys.exit(1)
                vcfSample = len(fd) - 9
                continue

            if ln[:3].upper() == "CHR" and not chrWarningGiven:
                chrWarningGiven = True

            vcfSiteLine += 1
            fd = ln.strip().split()
            if len(fd) != nField:
                print(f"Line [ {lineNo + 1} ] does not have correct column number, exiting!", file=logger)
                print(f"Current line has {len(fd)} columns.", file=logger)
                print(f"First 50 characters in the current line content [ {ln.strip()[:50]} ].", file=logger)
                sys.exit(1)

            chrom, pos, rsId, ref, alt, qual, filt, info, format = fd[:9]
            site = f"{chrom}:{pos}"

            # Check 1: Non-SNP
            if len(ref) != 1 or len(alt) != 1:
                fNonSnp.write(f"{' '.join(fd[:5])}\n")
                nNonSnp += 1
                if fExclude and rsId != '.':
                    fExclude.write(f"{rsId}\n")
                continue
            if ref not in ACGT or alt not in ACGTM:
                fNonSnp.write(f"{' '.join(fd[:5])}\n")
                nNonSnp += 1
                if fExclude and rsId != '.':
                    fExclude.write(f"{rsId}\n")
                continue

            # Check 2: Duplicates
            if site in snpSite:
                print(f"Duplicated site [ {site} ]", file=logger)
                fDup.write(f"DuplicatedSite\t{site}\n")
                nDupSite += 1
                if fExclude and rsId != '.':
                    fExclude.write(f"{rsId}\n")
                continue
            else:
                snpSite.add(site)

            # Check 3: Ascending order (just a warning, not excluded)
            if prevChrom != chrom:
                prevChrom, prevPos = chrom, int(pos)
            elif prevPos > int(pos):
                actionItem(logger)
                print(f"At line [ {lineNo + 1} ], genomic position {chrom}:{pos} is before previous position {prevChrom}:{prevPos}", file=logger)
                continue

            # Check 4: Reference mismatch
            trueRef = gs.getBase(chrom, int(pos) - 1)
            if trueRef is None:
                fRef.write(f"FailedGetBase\t{site}\n")
                nRef += 1
                if fExclude and rsId != '.':
                    fExclude.write(f"{rsId}\n")
                continue
            if ref != trueRef.decode('utf-8'):
                fRef.write(f"MismatchRefBase\t{site}:{trueRef.decode('utf-8')}-{ref}/{alt}\n")
                nRef += 1
                if fExclude and rsId != '.':
                    fExclude.write(f"{rsId}\n")
                continue

            # Check 5: Genotype format
            try:
                gtIndex = [idx for idx, i in enumerate(format.split(':')) if i == 'GT'][0]
            except IndexError:
                print(f"Line [ {lineNo} ] does not have GT defined in the FORMAT field", file=logger)
                continue

            try:
                genos = [i.split(':')[gtIndex] for i in fd[9:]]
            except IndexError:
                fGeno.write(f"IndividualMissingGTField\tLine:{lineNo + 1}\n")
                nGeno += 1
                if fExclude and rsId != '.':
                    fExclude.write(f"{rsId}\n")
                continue
            if not all([checkGTformat(g) for g in genos]):
                fGeno.write(f"IndividualHasInvalidGT\tLine:{lineNo + 1}\n")
                nGeno += 1
                if fExclude and rsId != '.':
                    fExclude.write(f"{rsId}\n")
                continue

            # Check 6: Monomorphic and AF > 0.5
            geno = [getGeno(g) for g in genos]
            ac = sum(g for g in geno if g > 0)
            nSample = sum(1 for g in geno if g >= 0)
            if ac == 0 or ac == 2 * nSample:
                fMono.write(f"{site}\t{ac}\t{nSample}\n")
                nMono += 1
                if fExclude and rsId != '.':
                    fExclude.write(f"{rsId}\n")

            if nSample > 0:
                af = 1.0 * ac / nSample / 2
            else:
                af = 0.0
            if af > 0.5:
                fAF.write(f"{site}\t{ref}\t{alt}\t{af}\n")
                nAF += 1
                if fExclude and rsId != '.':
                    fExclude.write(f"{rsId}\n")

    except SystemExit:
        sys.exit(1)
    except KeyboardInterrupt:
        print(f"VCF checking stopped at line [ {lineNo + 1} ]", file=logger)
        print(f" [ {ln[:50]} ... ]", file=logger)
        sys.exit(1)
    except Exception as e:
        print(f"VCF checking failed at line [ {lineNo + 1} ]", file=logger)
        print(f" [ {ln[:50]} ... ]", file=logger)
        print(f"Python exceptions occurred [ {e} ]!", file=logger)
        print("Please report the above to zhanxw@gmail.com", file=logger)
        raise
        sys.exit(1)

    # Final report
    if chrWarningGiven:
        print("---------------     WARNING     ---------------", file=logger)
        print("Detected that chromosome names have 'chr' prefix...", file=logger)
        print("Please consider using: '(grep ^\"#\" $your_old_vcf; grep -v ^\"#\" $your_old_vcf | sed 's:^chr::ig' | sort -k1,1n -k2,2n) | bgzip -c > $your_vcf_file'", file=logger)

    print("---------------     REPORT     ---------------", file=logger)
    print(f"Total [ {lineNo + 1} ] lines processed", file=logger)
    print(f"Examine [ {vcfHeaderLine} ] VCF header lines, [ {vcfSiteLine} ] variant sites, [ {vcfSample} ] samples", file=logger)
    print(f"[ {nDupSite} ] duplicated sites", file=logger)
    print(f"[ {nNonSnp} ] NonSNP sites outputted to [ {outPrefix}.check.nonSnp ]", file=logger)
    print(f"[ {nRef} ] Inconsistent reference sites outputted to [ {outPrefix}.check.ref ]", file=logger)
    print(f"[ {nGeno} ] Variant sites with invalid genotypes outputted to [ {outPrefix}.check.geno ]", file=logger)
    print(f"[ {nAF} ] Alternative allele frequency > 0.5 sites outputted to [ {outPrefix}.check.af ]", file=logger)
    print(f"[ {nMono} ] Monomorphic sites outputted to [ {outPrefix}.check.mono ]", file=logger)
    if excludeFile:
        print(f"All failing SNPs (rsIDs) written to [ {excludeFile} ]", file=logger)

    # Close files
    fDup.close()
    fRef.close()
    fNonSnp.close()
    fMono.close()
    fGeno.close()
    fAF.close()
    if fExclude:
        fExclude.close()

    actionItem(logger)
    if nDupSite > 0 or nRef > 0 or nGeno > 0:
        if nDupSite > 0:
            print("* Remove duplicated sites and rerun checkVCF.py", file=logger)
        if nRef > 0:
            print(f"* Read {outPrefix}.check.ref, ensure forward strand for autosomal sites", file=logger)
        if nGeno > 0:
            print(f"* Open {outPrefix}.check.geno, verify genotypes at listed lines", file=logger)
    else:
        print("* No errors found by checkVCF.py, thank you for cleaning the VCF file.", file=logger)
    print(f"* Upload these files to the FTP server for review: {outPrefix}.check.log {outPrefix}.check.dup {outPrefix}.check.nonSnp {outPrefix}.check.ref {outPrefix}.check.geno {outPrefix}.check.af {outPrefix}.check.mono" + (f" {excludeFile}" if excludeFile else ""), file=logger)
