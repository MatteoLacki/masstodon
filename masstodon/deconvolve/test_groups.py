from masstodon.data.substanceP_wh15_wv400 import mz, intensity, fasta, modifications, z
from masstodon.masstodon import masstodon_single


def test_groups():
    threshold = .025
    m, t = masstodon_single(mz, intensity, fasta, z,
                         modifications      = modifications,
                         orbitrap           = False,
                         threshold          = threshold,
                         isotopic_coverage  = .999,
                         min_prob           = .8, 
                         std_cnt            = 3,
                         get_timings        = True,
                         deconvolution_method = "nnls")

    assert all(m.imperator.groups.max_mz - m.imperator.groups.min_mz <= 
               threshold*2.00001)


def too_much_of_z4():
    threshold = .025
    m, t = masstodon_single(mz, intensity, fasta, z,
                         modifications      = modifications,
                         orbitrap           = False,
                         threshold          = threshold,
                         isotopic_coverage  = .999,
                         min_prob           = .8, 
                         std_cnt            = 3,
                         get_timings        = True,
                         deconvolution_method = "nnls")
    for M in m.ome.observables():
        if M.name == 'z4':
            z4_intensity = M.intensity
            break
    assert z4_intensity < 40000
