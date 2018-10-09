from masstodon.read.npy import spectrum_from_npy
from masstodon.plot.spectrum import plot_spectrum
import matplotlib.pyplot as plt

S1 = "/Users/matteo/Downloads/spec_comparison/S1"
S2 = "/Users/matteo/Downloads/spec_comparison/S2"
S3 = "/Users/matteo/Downloads/spec_comparison/S3"

S1 = spectrum_from_npy(S1)
S2 = spectrum_from_npy(S2)
S3 = spectrum_from_npy(S3)

len(S1[0])
len(S2[0])
len(S3[0])

# plot_spectrum(*S1, show=True)


plt.subplot(3, 1, 1)
plot_spectrum(*S1, show=False)
plt.subplot(3, 1, 2)
plot_spectrum(*S2, show=False)
plt.subplot(3, 1, 3)
plot_spectrum(*S3, show=False)
plt.show()
			
from masstodon.read.txt import spectrum_from_txt

s1 = spectrum_from_txt("/Users/matteo/Downloads/spec_comparison/S1/spec_1.txt")
s2 = spectrum_from_txt("/Users/matteo/Downloads/spec_comparison/S2/spec_2.txt")
s3 = spectrum_from_txt("/Users/matteo/Downloads/spec_comparison/S3/spec_3.txt")

plt.subplot(2, 1, 1)
plot_spectrum(*S1, show=False)
plt.subplot(2, 1, 2)
plot_spectrum(*s1, show=False)
plt.show()


plt.subplot(2, 1, 1)
plot_spectrum(*S2, show=False)
plt.subplot(2, 1, 2)
plot_spectrum(*s2, show=False)
plt.show()

plt.subplot(2, 1, 1)
plot_spectrum(*S3, show=False)
plt.subplot(2, 1, 2)
plot_spectrum(*s3, show=False)
plt.show()


# ./printspectrum -sn 1 -raw ~/Projects/masstodon/data/PXD001845/raw_files/20141202_AMB_Bora_10x_40MeOH_1FA_OT_120k_10uscans_795_ETD_4ms_22precZ.raw > ~/Downloads/spec_comparison/S1/spec_1.txt
