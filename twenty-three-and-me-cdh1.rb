#!/usr/bin/env ruby
#------------------------------------------------------------------------------#
#                        twenty-three-and-me-cdh1.rb                           #
#------------------------------------------------------------------------------#

# Check a raw data file from 23 and me for SNPs in the CDH1 (E-cadherin) gene




#---------------------------- Require libraries -------------------------------#

require 'optparse'
require 'set'




#------------------ Initialize command line option struct ---------------------#

Options = Struct.new(
  :bed_file_path,
  :fasta_file_path,
  :vcf_file_path,
  :raw_file_path,
  :annotation
)




#---------------------------- Class definitions -------------------------------#

# A parser for command line arguments.
class Parser
  def self.parse(options)
    args = Options.new('world')

    opt_parser = OptionParser.new do |parser|
      parser.banner = 'Usage: twenty-three-and-me-cdh1.rb -b BED -f FASTA ' \
        '-v VCF -r RAW [-a ANNOTATION]'

      parser.on(
        '-b',
        '--bed BED',
        'Path to BED file with gene annotations'
      ) do |bed_file_path|
        args.bed_file_path = bed_file_path
      end
      
      parser.on(
        '-f',
        '--fasta FASTA',
        'Path to FASTA file containing reference CDH1 region'
      ) do |fasta_file_path|
        args.fasta_file_path = fasta_file_path
      end
      
      parser.on(
        '-v',
        '--vcf VCF',
        'Path to VCF file containing info on CDH1 region'
      ) do |vcf_file_path|
        args.vcf_file_path = vcf_file_path
      end
      
      parser.on(
        '-r',
        '--23-and-me-raw RAW',
        'Path to 23 and me raw data file'
      ) do |raw_file_path|
        args.raw_file_path = raw_file_path
      end
      
      parser.on(
        '-a',
        '--annotation ANNOTATION',
        'Limit output to variants with the provided annotation'
      ) do |annotation|
        args.annotation = annotation
      end

      parser.on_tail('-h', '--help', 'Show this message') do
        puts parser
        exit
      end
    end

    opt_parser.parse!(options)
    args
  end
end

# A genetic variant or SNP.
class Variant
  def initialize(chromosome, position, id, genotype)
    @chromosome = chromosome
    @position = position
    @id = id
    @genotype = genotype
    @annotations = Set.new
  end
  
  def check_annotation(annotation, interval)
    start, stop = interval
    if @position >= start && @position <= stop
      @annotations.add(annotation)
    end
  end
  
  def add_annotation(annotation)
    @annotations.add(annotation)
  end
  
  def check_reference_allele(region_start, region_sequence)
    @reference_allele = region_sequence[@position - region_start - 1]
  end
  
  def chromosome
    @chromosome
  end
  
  def position
    @position
  end
  
  def id
    @id
  end
  
  def genotype
    @genotype
  end
  
  def annotations
    @annotations
  end
  
  def region_start
    @region_start
  end
  
  def reference_allele
    @reference_allele
  end
end



#---------------------------- method definitions ------------------------------#

# Loads the gene annotations.
#
# file_path - String giving the path to the gene annotations file (BED format)
#             on disk.
#
# Returns an array containing the coordinates of each exon and intron.
def load_gene_annotations(file_path)
  annotations = { 'exon' => [], 'intron' => [] }
  
  File.open(file_path).each do |line|
    chromosome, start, stop, annotation = line.split
    annotations[annotation] = [start.to_i, stop.to_i]
    REGION_TYPES.each do |region_type|
      if annotation == nil
        next
      elsif annotation.start_with?(region_type)
        annotations[region_type] << annotations[annotation]
      end
    end
  end
  
  annotations
end

# Loads sequence data from the FASTA file.
#
# file_path - String giving the path to the FASTA file containing the gene
#             sequence.
#
# Returns an array representation of the FASTA file.
def load_fasta(file_path)
  fasta = {}
  region = ''
  
  File.open(file_path).each do |line|
    if line.start_with?('>')
      region = line[1..-2]
    else
      fasta[region] = line.rstrip
    end
  end
  
  fasta
end

# Loads variant data from the VCF file.
#
# file_path - String giving the path to the VCF file containing 1000 Genomes
#             variant information (e.g. allele frequencies).
#
# Returns an array representation of the VCF file.
def load_vcf(file_path)
  vcf = {}
  
  File.open(file_path).each do |line|
    unless line.start_with?('#')
      chrom, pos, id, ref, alt, _, _, info = line.split
      vcf[id] = { chrom: chrom, pos: pos.to_i, ref: ref, alt: alt, info: {} }
      info.split(';').each do |entry|
        key, value = entry.split('=')
        vcf[id][:info][key] = value
      end
    end
  end
  
  vcf
end

# Scans the raw data file for variants in the CDH1 region and processes the
# variants.
#
# file_path   - String giving the path to the 23andMe raw data file
# annotations - Array containing gene annotations
# fasta       - Array containing gene sequence
#
# Returns an array of variants falling within the CDH1 region and their
# genotypes per the raw data file.
def scan_raw_data_file(file_path, annotations, fasta)
  variants = []
  
  File.open(file_path).each do |line|
    id, chromosome, position, genotype = line.split
    position = position.to_i
    if [
      chromosome == '16',
      position >= annotations['promoter'][0],
      position <= annotations['3UTR'][1]
    ].all?
      variant = Variant.new(chromosome, position, id, genotype)
      
      annotations.each do |key, value|
        if REGION_TYPES.include?(key)
          annotations[key].each do |interval|
            variant.check_annotation(key, interval)
          end
        else
          variant.check_annotation(key, value)
          if variant.annotations.include?(key)
            start, stop = value
            variant.check_reference_allele(start, fasta[key])
          end
        end
      end
      
      variants << variant
    end
  end
  
  variants
end

# Checks which variants have uncommon (possibly delterious) alleles, and
# updates their variant annotations accordingly.
#
# variants - Array of variants extracted from 23andMe raw data file.
# vcf      - Array containing 1000 genomes variant information.
def check_for_uncommon_alleles(variants, vcf)
  variants.each do |variant|
    if variant.genotype == '--' || vcf[variant.id] == nil
      next
    elsif variant.genotype[0] != variant.genotype[1]
      variant.add_annotation('uncommon_allele')
    elsif [
      variant.genotype[0] != variant.reference_allele.upcase,
      variant.reference_allele.upcase == vcf[variant.id][:ref],
      vcf[variant.id][:info]['EUR_AF'].to_f < 0.5
    ].all?
      variant.add_annotation('uncommon_allele')
    elsif [
      variant.genotype[0] != variant.reference_allele.upcase,
      variant.reference_allele.upcase == vcf[variant.id][:alt],
      vcf[variant.id][:info]['EUR_AF'].to_f > 0.5
    ].all?
      variant.add_annotation('uncommon_allele')
    elsif [
      variant.genotype[0] == variant.reference_allele.upcase,
      variant.reference_allele.upcase == vcf[variant.id][:ref],
      vcf[variant.id][:info]['EUR_AF'].to_f <= 0.5
    ].all?
      variant.add_annotation('common_allele')
    elsif [
      variant.genotype[0] == variant.reference_allele.upcase,
      variant.reference_allele.upcase == vcf[variant.id][:alt],
      vcf[variant.id][:info]['EUR_AF'].to_f >= 0.5
    ].all?
      variant.add_annotation('common_allele')
    elsif [
      variant.genotype[0] != variant.reference_allele.upcase,
      variant.reference_allele.upcase == vcf[variant.id][:alt],
      vcf[variant.id][:info]['EUR_AF'].to_f <= 0.5
    ].all?
      variant.add_annotation('common_allele')
    elsif [
      variant.genotype[0] != variant.reference_allele.upcase,
      variant.reference_allele.upcase == vcf[variant.id][:ref],
      vcf[variant.id][:info]['EUR_AF'].to_f >= 0.5
    ].all?
      variant.add_annotation('common_allele')
    end
  end
end

# Checks if any of the variants mutate the protein sequence, and annotates them
# accordingly.
#
# variants    - Array of variants extracted from the 23andMe raw data file.
# annotations - Array containing gene annotations
# fasta       - Array containing gene sequence
def check_for_protein_mutations(variants, annotations, fasta)
  variants.each do |variant|
    unless (variant.annotations.include? 'exon') &&
             (
               variant.genotype[0] != variant.reference_allele.upcase ||
                 variant.genotype[1] != variant.reference_allele.upcase
             )
      next
    end
      
    if variant.genotype[0] != variant.reference_allele.upcase
      alternate_allele = variant.genotype[0]
    elsif variant.genotype[1] != variant.reference_allele.upcase
      alternate_allele = variant.genotype[1]
    end
    
    exon=''
    variant.annotations.each do |annotation|
      if annotation.start_with?('exon_')
        exon = annotation
      end
    end
    
    exon_start = annotations[exon][0]
    sequence = fasta[exon]
    sequence_position = variant.position - exon_start - 1
    codon_position = sequence_position % 3
    codon = sequence.scan(/.{3}/)[(sequence_position - codon_position) / 3]
              .upcase
    variant_codon = String.new(codon)
    variant_codon[codon_position] = alternate_allele
    
    if GENETIC_CODE[codon] == GENETIC_CODE[variant_codon]
      next
    elsif GENETIC_CODE[variant_codon] == 'stop'
      variant.add_annotation('nonsense')
    else
      variant.add_annotation('missense')
    end
  end
end

# Main method
def main(options)
  annotations = load_gene_annotations(options.bed_file_path)
  fasta = load_fasta(options.fasta_file_path)
  vcf = load_vcf(options.vcf_file_path)
  variants = scan_raw_data_file(options.raw_file_path, annotations, fasta)
  check_for_uncommon_alleles(variants, vcf)
  check_for_protein_mutations(variants, annotations, fasta)
  
  puts [
    '#rsid',
    'chromosome',
    'position',
    'genotype',
    'reference_allele',
    'european_alt_allele_frequency',
    'annotations'
  ].join("\t")
  
  variants.each do |variant|
    unless (variant.annotations.include?(options.annotation)) ||
             (options.annotation == nil)
      next
    end
      
    puts [
      variant.id,
      variant.chromosome,
      variant.position.to_s,
      variant.genotype,
      variant.reference_allele.upcase,
      vcf[variant.id] ? vcf[variant.id][:info]['EUR_AF'] : 'nil',
      variant.annotations.to_a.join(',')
    ].join("\t")
  end
end

# Parses the command-line options
def parse_options()
  Parser.parse ARGV
end




#--------------------------------- constants ----------------------------------#

REGION_TYPES = %w(exon intron)

GENETIC_CODE = {
  'AAA' => 'Lys',
  'AAC' => 'Asn',
  'AAG' => 'Lys',
  'AAT' => 'Asn',
  'ACA' => 'Thr',
  'ACC' => 'Thr',
  'ACG' => 'Thr',
  'ACT' => 'Thr',
  'AGA' => 'Arg',
  'AGC' => 'Ser',
  'AGG' => 'Arg',
  'AGT' => 'Ser',
  'ATA' => 'Ile',
  'ATC' => 'Ile',
  'ATG' => 'Met',
  'ATT' => 'Ile',
  'CAA' => 'Gln',
  'CAC' => 'His',
  'CAG' => 'Gln',
  'CAT' => 'His',
  'CCA' => 'Pro',
  'CCC' => 'Pro',
  'CCG' => 'Pro',
  'CCT' => 'Pro',
  'CGA' => 'Arg',
  'CGC' => 'Arg',
  'CGG' => 'Arg',
  'CGT' => 'Arg',
  'CTA' => 'Leu',
  'CTC' => 'Leu',
  'CTG' => 'Leu',
  'CTT' => 'Leu',
  'GAA' => 'Glu',
  'GAC' => 'Asp',
  'GAG' => 'Glu',
  'GAT' => 'Asp',
  'GCA' => 'Ala',
  'GCC' => 'Ala',
  'GCG' => 'Ala',
  'GCT' => 'Ala',
  'GGA' => 'Gly',
  'GGC' => 'Gly',
  'GGG' => 'Gly',
  'GGT' => 'Gly',
  'GTA' => 'Val',
  'GTC' => 'Val',
  'GTG' => 'Val',
  'GTT' => 'Val',
  'TAA' => 'stop',
  'TAC' => 'Tyr',
  'TAG' => 'stop',
  'TAT' => 'Tyr',
  'TCA' => 'Ser',
  'TCC' => 'Ser',
  'TCG' => 'Ser',
  'TCT' => 'Ser',
  'TGA' => 'stop',
  'TGC' => 'Cys',
  'TGG' => 'Trp',
  'TGT' => 'Cys',
  'TTA' => 'Leu',
  'TTC' => 'Phe',
  'TTG' => 'Leu',
  'TTT' => 'Phe'
}




#---------------------------------- execute -----------------------------------#

options = parse_options
main(options)
