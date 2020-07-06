##############################################
# pyplot
pyplot()
import PyCall
import PyPlot: matplotlib
using LaTeXStrings

# font
PyCall.PyDict(matplotlib["rcParams"])["font.family"] = "serif:bold"
# PyCall.PyDict(matplotlib["rcParams"])["font.weight"] = "bold"
PyCall.PyDict(matplotlib["rcParams"])["font.serif"] = "Helvetica" #"STIX"
PyCall.PyDict(matplotlib["rcParams"])["mathtext.fontset"] = "dejavusans"#"Helvetica"
PyCall.PyDict(matplotlib["rcParams"])["text.usetex"] = true
# PyCall.PyDict(matplotlib["rcParams"])["text.latex.unicode"] = true
# PyCall.PyDict(matplotlib["rcParams"])["text.latex.preamble"] = ["\\usepackage{siunitx}"];
