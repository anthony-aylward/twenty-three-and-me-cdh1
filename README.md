# twenty-three-and-me-cdh1
A ruby script to check a 23andMe raw data file for variants in the CDH1 (E-cadherin) gene

```
Usage: ruby twenty-three-and-me-cdh1.rb -b BED -f FASTA -v VCF -r RAW [-a ANNOTATION]
    -b, --bed BED                    Path to BED file with gene annotations
    -f, --fasta FASTA                Path to FASTA file containing reference CDH1 region
    -v, --vcf VCF                    Path to VCF file containing info on CDH1 region
    -r, --23-and-me-raw RAW          Path to 23 and me raw data file
    -a, --annotation ANNOTATION      Limit output to variants with the provided annotation
    -h, --help                       Show this message
```
