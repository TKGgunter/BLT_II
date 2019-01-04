import ROOT, ctypes

f = ROOT.TFile("egamma_tightSF.root")
graph = f.Get("EGamma_SF2D")

lines = ""
for j in range(1,8):
  for i in range(1,10):
    #print i, j, graph.GetBinContent(i,j)
    lines += str(graph.GetBinContent(i,j)) + " & "
  lines += "\/\ \n"

print lines

