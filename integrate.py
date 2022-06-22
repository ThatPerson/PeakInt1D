import argparse as ap
import nmrglue as ng
import matplotlib.pyplot as plt
from lmfit.models import VoigtModel, QuadraticModel, GaussianModel, LorentzianModel
import numpy as np

parser = ap.ArgumentParser()
parser.add_argument('-spectra', type=str, default=None, help='Spectra directory.')
parser.add_argument('-initial', type=str, default=None, required=True, help="Initial peak guesses")
parser.add_argument('-output', type=str, default="output.txt", help="Output text file.")
parser.add_argument('-vclist', required=True, type=str, help="vd/vplist file")
parser.add_argument('-vc', required=False, default=None, type=float, help="vc list cycle time")
parser.add_argument('-figs', required=False, default=None, type=str, help="output fig dir")
args = parser.parse_args()

def add_peak(prefix, center, amplitude=5e9, sigma=0.05, allow_neg=False):
	peak = LorentzianModel(prefix=prefix)
	pars = peak.make_params()
	pars[prefix + 'center'].set(center)
	if (allow_neg == False):
		pars[prefix + 'amplitude'].set(amplitude, min=0)
	else:
		pars[prefix + 'amplitude'].set(amplitude)
	pars[prefix + 'sigma'].set(sigma, min=0, max=0.7)
	#pars[prefix + 'gamma'].set(sigma, vary=True, min=0)
	return peak, pars

def rough_peaks(st):
	k = [float(p) for p in st.split(",")]
	model = QuadraticModel(prefix='bkg_')
	params = model.make_params(a=0, b=0, c=0)
	for i, cen in enumerate(k):
		allow_neg = False
		if (i == 0):	
			allow_neg = True
		peak, pars = add_peak('vz%d_' % (i+1), cen, allow_neg=allow_neg)
		model = model + peak
		params.update(pars)
	return model, params

def read_spectra(fn, vclist, vcstep):
	dic, data = ng.bruker.read_pdata(fn)
	F1P = dic['procs']['F1P']
	F2P = dic['procs']['F2P']
	k = np.shape(data)
	xscale = np.linspace(F1P, F2P, k[1])
	dat = []
	vcl = np.loadtxt(vclist) * vcstep
	if (len(vcl) >= k[0]):
		vcl = vcl[:k[0]]
	else:
		print("Incorrect length VC list")
		exit(-1)
		
	for i in range(0, k[0]):
		dat.append((vcl[i], data[i, :]))
	
	return dat, xscale

dat, xscale = read_spectra(args.spectra, args.vclist, args.vc)
mdl, prms = rough_peaks(args.initial)
init = mdl.eval(prms, x=xscale)
ylims = None
with open(args.output, "w") as f:
	header = True
	for time, spec in dat:
		result = mdl.fit(spec, prms, x=xscale)
		comps = result.eval_components()
		#result.conf_interval()
		print(result.fit_report(min_correl=0.5))
		if (args.figs is not None):
			plt.plot(xscale, spec, label='data')
			plt.plot(xscale, result.best_fit, label='best fit')
			for name, comp in comps.items():
				plt.plot(xscale, comp, '--', label=name)
				print("%s: %f" % (name, np.sum(comp) / np.sum(result.best_fit)))
			plt.legend(loc='upper right')
			plt.xlim([17.5, -5])
			if (ylims is None):
				ylims = plt.gca().get_ylim()
			else:
				plt.gca().set_ylim(ylims)
			plt.savefig("%s/%0.3f.png" % (args.figs, time))
			plt.close()
		if (header):
			f.write("# Time")
			for name, comp in comps.items():
				f.write(", %s" % (name))
			f.write("\n")
			header = False
		f.write("%0.6f" % (time))
		print(result.params)
		print(result.ci_out)
		for i in range(0, len(args.initial.split(","))):
			val = result.params['vz%d_amplitude' % (i+1)].value
			stderr = result.params['vz%d_amplitude' % (i+1)].stderr
			if (stderr is None):
				stderr = -1
			f.write(", %f, %f" % (val, stderr))
		f.write("\n")
		
		for prm in prms:
			if ("center" in prm):
				prms[prm].set(value=result.params[prm].value, vary=False)
#				prms.set(prm, vary=False)
				print(prm)
		
