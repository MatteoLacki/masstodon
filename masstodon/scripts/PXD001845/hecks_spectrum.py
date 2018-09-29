from os.path import join as pjoin

from masstodon.read.npy  import spectrum_from_npy



def get_hecks_spectrum(all_data, folder, scan):
	data_path = pjoin(all_data, folder, str(scan))
	return spectrum_from_npy(data_path)