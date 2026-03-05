# =============================================================================
# atlas.R
# Harvard-Oxford Atlas (HOA) — 110 region MNI coordinates and lobe assignments
# Coordinates sourced from braingraph library via paper appendix (test.tex).
# Used for anatomical node placement in brain network visualizations.
# =============================================================================

HOA_ATLAS <- data.frame(
  name = c(
    "Thal.L","Thal.R","Caud.L","Caud.R","Put.L","Put.R","Pall.L","Pall.R",
    "Hip.L","Hip.R","Amy.L","Amy.R","Accbns.L","Accbns.R",
    "FP.L","FP.R","INS.L","INS.R","F1.L","F1.R","F2.L","F2.R",
    "F3t.L","F3t.R","F3o.L","F3o.R","PRG.L","PRG.R",
    "TP.L","TP.R","T1a.L","T1a.R","T1p.L","T1p.R",
    "T2a.L","T2a.R","T2p.L","T2p.R","TO2.L","TO2.R",
    "T3a.L","T3a.R","T3p.L","T3p.R","TO3.L","TO3.R",
    "POG.L","POG.R","SPL.L","SPL.R","SGa.L","SGa.R","SGp.L","SGp.R",
    "AG.L","AG.R","OLs.L","OLs.R","OLi.L","OLi.R",
    "CALC.L","CALC.R","FMC.L","FMC.R","SMC.L","SMC.R",
    "SC.L","SC.R","PAC.L","PAC.R","CGa.L","CGa.R","CGp.L","CGp.R",
    "PCN.L","PCN.R","CN.L","CN.R","FOC.L","FOC.R",
    "PHa.L","PHa.R","PHp.L","PHp.R","LG.L","LG.R",
    "TFa.L","TFa.R","TFp.L","TFp.R","TOF.L","TOF.R","OF.L","OF.R",
    "FO.L","FO.R","CO.L","CO.R","PO.L","PO.R",
    "PP.L","PP.R","H.L","H.R","PT.L","PT.R",
    "SCLC.L","SCLC.R","OP.L","OP.R"
  ),
  x = c(
    -9.99,10.92,-12.84,13.5,-24.79,25.48,-19.13,19.75,
    -25.18,26.5,-22.86,22.77,-9.34,9.21,
    -25.01,26.44,-36.42,37.5,-14.56,15.09,-38.15,39.19,
    -49.76,51.84,-50.69,52.49,-34.28,35.08,
    -40.44,40.99,-56.0,57.22,-62.37,61.35,
    -57.8,57.86,-60.95,60.97,-57.4,58.32,
    -47.97,46.31,-53.46,53.8,-51.81,54.19,
    -38.53,37.25,-29.28,29.08,-57.0,58.23,-54.86,55.23,
    -50.45,52.15,-32.05,33.03,-45.21,45.47,
    -10.47,11.93,-5.38,5.51,-5.79,6.39,
    -5.7,5.66,-6.82,7.07,-5.15,5.75,-6.3,6.96,
    -8.16,9.33,-8.66,9.38,-29.69,29.31,
    -21.68,22.59,-22.15,22.95,-12.57,13.93,
    -32.3,30.87,-36.02,36.53,-33.32,35.02,-26.33,27.24,
    -39.82,41.14,-48.03,49.47,-48.4,48.86,
    -46.77,48.11,-45.22,46.04,-52.64,54.84,
    -12.29,9.09,-17.23,18.09
  ),
  y = c(
    -19.16,-18.5,8.92,9.66,0.54,2.03,-5.17,-3.87,
    -23.25,-20.99,-5.18,-3.69,11.14,11.46,
    52.77,52.0,1.01,2.65,17.96,17.55,18.31,18.3,
    28.6,27.82,14.63,15.48,-11.71,-10.59,
    11.07,12.93,-3.79,-0.98,-29.14,-23.87,
    -4.41,-1.74,-27.39,-22.35,-52.7,-49.3,
    -5.1,-2.16,-28.21,-23.36,-53.45,-49.71,
    -27.78,-26.58,-49.4,-47.79,-32.5,-27.26,-46.04,-40.29,
    -55.74,-51.69,-72.77,-71.06,-75.63,-74.11,
    -74.92,-73.63,43.84,43.41,-2.67,-2.85,
    20.6,20.42,36.58,36.37,18.19,19.31,-38.56,-35.8,
    -60.06,-58.48,-80.04,-78.23,23.81,23.43,
    -9.28,-8.04,-32.31,-30.25,-65.51,-62.73,
    -4.53,-2.55,-29.45,-23.81,-53.65,-49.88,-76.86,-75.48,
    18.29,18.82,-8.28,-5.66,-31.53,-27.69,
    -5.34,-3.56,-20.04,-17.36,-29.69,-25.33,
    -68.86,-74.03,-96.36,-95.14
  ),
  z = c(
    16.28,16.6,19.71,20.87,10.42,10.35,8.53,8.56,
    -3.92,-4.07,-7.49,-7.91,2.89,3.57,
    17.75,18.6,10.16,9.83,66.57,67.52,52.01,53.02,
    18.59,17.72,25.22,26.37,59.18,59.79,
    -19.78,-19.31,1.86,-0.41,13.86,11.5,
    -12.05,-14.52,-0.91,-2.18,10.87,11.53,
    -29.12,-31.18,-16.0,-18.1,-6.68,-6.86,
    61.5,62.93,67.63,68.92,46.94,48.18,43.58,43.9,
    39.3,42.16,47.99,49.0,8.06,8.49,
    18.19,18.36,-7.87,-8.24,66.3,67.64,
    -5.68,-5.9,30.93,32.84,34.6,34.15,38.79,40.04,
    47.25,48.1,37.62,37.98,-6.49,-6.21,
    -20.7,-20.63,-7.13,-6.98,4.55,5.03,
    -31.6,-32.28,-15.04,-18.05,-5.95,-6.56,-3.45,-2.29,
    14.62,14.75,21.64,21.1,30.3,31.65,
    2.46,2.9,17.31,16.89,20.8,22.39,
    25.57,24.47,17.23,18.35
  ),
  lobe = c(
    "SCGM","SCGM","SCGM","SCGM","SCGM","SCGM","SCGM","SCGM",
    "SCGM","SCGM","SCGM","SCGM","SCGM","SCGM",
    "Frontal","Frontal","Insula","Insula","Frontal","Frontal","Frontal","Frontal",
    "Frontal","Frontal","Frontal","Frontal","Frontal","Frontal",
    "Temporal","Temporal","Temporal","Temporal","Temporal","Temporal",
    "Temporal","Temporal","Temporal","Temporal","Temporal","Temporal",
    "Temporal","Temporal","Temporal","Temporal","Temporal","Temporal",
    "Parietal","Parietal","Parietal","Parietal","Parietal","Parietal","Parietal","Parietal",
    "Parietal","Parietal","Occipital","Occipital","Occipital","Occipital",
    "Occipital","Occipital","Frontal","Frontal","Frontal","Frontal",
    "Frontal","Frontal","Cingulate","Cingulate","Cingulate","Cingulate","Cingulate","Cingulate",
    "Parietal","Parietal","Occipital","Occipital","Frontal","Frontal",
    "Temporal","Temporal","Temporal","Temporal","Occipital","Occipital",
    "Temporal","Temporal","Temporal","Temporal","Occipital","Occipital","Occipital","Occipital",
    "Frontal","Frontal","Frontal","Frontal","Parietal","Parietal",
    "Temporal","Temporal","Temporal","Temporal","Temporal","Temporal",
    "Occipital","Occipital","Occipital","Occipital"
  ),
  stringsAsFactors = FALSE
)

# Lobe colour palette for node colouring
LOBE_COLOURS <- c(
  Frontal   = "#E74C3C",
  Temporal  = "#3498DB",
  Parietal  = "#2ECC71",
  Occipital = "#F39C12",
  Cingulate = "#9B59B6",
  Insula    = "#1ABC9C",
  SCGM      = "#95A5A6"
)
