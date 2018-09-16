%load_ext autoreload
%autoreload 2
%load_ext line_profiler

from time import time

from masstodon.masstodon import Masstodon, masstodon_batch, masstodon_single, masstodon_load
from masstodon.ome.ome   import Ome
from masstodon.readers.from_npy import spectrum_from_npy

data_path     = '/Users/matteo/Projects/masstodon/data/PXD001845/numpy_files/20141202_AMB_pBora_PLK_10x_40MeOH_1FA_OT_120k_10uscans_928_ETciD_8ms_15SA_19precZ/1'
mz, intensity = spectrum_from_npy(data_path)

fasta             = "GAASMMGDVKESKMQITPETPGRIPVLNPFESPSDYSNLHEQTLASPSVFKSTKLPTPGKFRWSIDQLAVINPVEIDPEDIHRQALYLSHSRIDKDVEDKRQKAIEEFFTKDVIVPSPWTDHEGKQLSQCHSSKCTNINSDSPVGKKLTIHSEKSD"
charge            = 24
min_prob          = .8
isotopic_coverage = .999
std_cnt           = 3

precursors = [{'fasta': fasta,
               'name' : 'Abudhabizine',
               'q'    : 24}]

# m = masstodon_single(mz, intensity)
m = masstodon_batch(mz, intensity, precursors, 
                    isotopic_coverage=isotopic_coverage,
                    min_prob=min_prob, 
                    std_cnt=std_cnt)
n = masstodon_load("dump")
o = masstodon_single(mz, intensity, fasta, charge, "Yani")
p = masstodon_single(mz, intensity, fasta, charge, "Yani",
                     deconvolution_graph_path="dump/deconvolution_graph.gpickle")




# precursors = [
#     dict(fasta = "AAAACCCKKK",
#          name  = "a-zine",
#          q     = 1),
#     dict(fasta = "VTRSQLM",
#          q     = 2,
#          name  = "b-zine"),
#     dict(fasta = "QQQQNNNKKKHIFFDDD",
#          q     = 2,
#          name  = "c-zine"),
# ]
# molecules = [
#     dict(name   = 'saliva',
#          formula= 'C100H202',
#          q      =  1),
#     dict(name   = 'shit',
#          formula= 'C102H205',
#          q      =  2),
#     dict(name   = 'lipstick',
#          formula= 'C100H202',
#          q      =  2),
#     dict(name   = "startrek's leather underwear",
#          formula= 'C100H202',
#          q      =  2)
# ]


m = masstodon_batch(mz, intensity, precursors, 
                    isotopic_coverage=isotopic_coverage,
                    min_prob=min_prob, 
                    std_cnt=std_cnt)




todon = Masstodon(mz        = mz, 
                  intensity = intensity, 
                  mz_digits = 4,
                  precursors= precursors,
                  # molecules = molecules,
                  # fasta     = fasta,
                  # charge    = charge,
                  # name             = "", 
                  # modifications    = {},
                  # fragments        = "cz",
                  # blocked_fragments= ['c0'],
                  # block_prolines   = True,
                  # distance_charges = 5.,
                  std_cnt          = 3,
                  isotopic_coverage= .999,
                  min_prob         = .7,
                  deconvolution_graph_path = '')

todon.set_spectrum()
todon.set_isotopic_calculator()
todon.set_ome()
# todon.ome.plot()


t0 = time()
ome = Ome(todon.iso_calc)
ome.make_graph(precursors)
t1 = time()

good_mols, good_subspectra = ome.filter_by_deviations(todon.subspectra, todon.std_cnt)





ome = Ome(todon.iso_calc)
ome.make_graph(precursors, molecules)
ome.plot()
list(ome.observables())
list(ome.sources())

molecules = [
    dict(name   = 'saliva',
         formula= 'C100H202',
         q      =  2),
    dict(name   = 'lipstick',
         formula= 'C100H202',
         q      =  2),
    dict(name   = "startrek's leather underwear",
         formula= 'C100H202',
         q      =  2)
]
ome = Ome(todon.iso_calc)
ome.make_graph([], molecules)
ome.plot()



# still doesn't work. 
# we need to be able to distinguish the source nodes.


for n in ome.G:
    print(n)
    print(n.__hash__())

