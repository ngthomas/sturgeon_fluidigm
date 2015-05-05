

# These are just some notes I am compiling while doing some literature surveys:


## Analysing polyploids via dominant markers:
* Falush et al 2007 (AFLP structure paper) the limitations of genotype calling
in polyploids mean that many conventional analysis methods are invalid for many organisms.
They use a dominant approach for population clustering
* Wang and Scribner.  Doing the same for parentage inference.
    * They note that a common approach for polyploids is to make dominant pseudo loci out of the genotypes
    * Their abstract makes it clear that they are focusing on microsats here. They don't really mention SNPs at all.  
    * _In the discussion they talk a little about inheritance modes in a nice way: In reality, most polyploid species
    may fall somewhere on a continuum between the two extreme inheritance modes. Autopolyploids result from genome
    duplication events and should in principle follow the nondisomical inheritance. However, accrual of mutational
    differences in microsatellite motif size and at priming sites (Estoup & Cornuet 1999) can lead to loss of
    amplification prod- ucts in different gene copies and to an apparently diso- mical inheritance of markers.
    Allopolyploids are due to hybridization events and thus should in principle fol- low the disomical inheritance.
    However, the hybridiza- tions are typically between closely related species with similar genomes, and thus,
    allopolyploids may still dis- play nondisomical inheritance to some extent. 
    * They have a long discussion of estimating selfing rates using the Rodzen transformation and point out that
    it would be tough to to it with the program inStruct with polyploid data.  You would have a hard time doing it
    with genotypes called the way we have as well.
* Rodzen JA, Famula TR, May B (2004) Estimation of parentage and relatedness in the polyploid white sturgeon (Acipenser transmontanus) using a dominant marker approach for duplicated microsatellite loci. _This paper is cited
by Wang and Scribner as one of the early formal application of the pseudodominance method, though others were used
earlier, too_  It will be worth citing and discussion this, since our method is quite different---_because we can't 
call individual alleles all that reliably_.


## Other software for analyzing polyploids
*  Puyvelde and friends:  ATETRA.  Basically trying to infer allele frequencies from partially observed heterozygotes.
Totally microsat based.


## SNP discovery and typing of polyploids
* Rob Ogden's paper has a good overview. He would be a good guy to send our paper to. In his paper, their goal was to find
  SNPs that could distinguish lots of different populations of sturgeons.
    * They had a family group in there to try to assess Mendelian inheritance. Although that was maybe
    more for their discovery phase).
    * A total of 675 candidate SNP markers showed alternate fixed alleles between Russian and Persian
    samples in the sequence data (Table S4.2, Supporting information). Of those carried forward for
    population testing, three of the five Kaspar genotyping assays generated reproduc- ible data with
    one showing clustering into the two putative species (Fig. 3); the remaining two displayed poor 
    segregation among Russian and Caspian samples. _Take-home message: It is hard to find good SNPs in polyploid sturgeon!_
    * (Note to self.  With our method for pop assignment, you can not worry about inheritance so much.)
    * _They note again that genotyping of tetraploids has been tough_: In tetraploid species, with intermediate heterozygote
    genotypes, the ability to genotype at all has been severely restricted, although recent analytical approaches have been
    designed to address this issue for SNP genotyping assays (e.g. Serang et al. 2012). 
* Seran and Garcia for typing SNPs from assays.
* Voorrips et al 2011: the tetraFit program.
    * Has a nice overview of how modern SNP assays work and why they are difficult to deal with for tetraploid data.
    Our main point is to say that it is there, but that it might be of limited utility when things can't be reliably
    scored using it.  (Show a two-cluster genotype and ask how tetraFit might try to call genotypes).
    * They do mention the difficulty of null alleles.  And, in our case, I think a lot of the "2-cluster" loci that
    are quite diagnostic have a null allele at very high frequency in one of the populations.  So, I will have to
    work in a discussion of null alleles.  They note that they can't really deal very well with null alleles.
* Gidskehaug L, Kent M, Hayes BJ, Lien S: Genotype calling and mapping of multisite variants using an Atlantic
salmon iSelect SNP array. Bioinformatics 2011, 27:303-310.  An interesting paper with a method for Illumina SNP chip
stuff for organisms with a variety of different 


## Generic SNP genotype calling by clustering
*[Takiotoh and friends](http://bioinformatics.oxfordjournals.org.oca.ucsc.edu/content/23/4/408.abstract) have worked
on generic clustering and automatic genotype calling from 2-D clusters.  But we will say that hand scoring for
96 SNPs is just fine for us.


## SNPs, Tetraploids and population structure
* Stich et al 2013: tetraploid potatos typed on Illumina Infinium chip.  46 individuals at 8303 SNPs scored _by hand_!
because they can't call genotypes automatically.  They used a 5-genotype system.  Then they did their population
structure clustering using principle coordinate analysis.  That seems typical.  I think it would be worthwhile to
develop a genotype-based structure analysis to see if the outliers (in our sturgeon data) in PCA look like a
separate cluster or not.