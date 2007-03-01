(TeX-add-style-hook "elementlibmanual"
 (lambda ()
    (LaTeX-add-labels
     "Truss1d"
     "Truss2d"
     "beam2dfig"
     "beam3dfig"
     "Planestress2dfig"
     "qplanstrssfig"
     "TrPlanestressfig"
     "qtrplanstressfig"
     "Quad1PlaneStrainfig"
     "TrplaneStrain"
     "Quad1ht"
     "Quad1htfig"
     "Tr1ht"
     "Tr1htfig"
     "Brick1ht"
     "Brick1htfig"
     "Tr1CBS"
     "Tr1CBSfig")
    (TeX-add-symbols
     '("param" 1)
     '("optelemparam" 2)
     '("elemparam" 2)
     '("elemkeyword" 1)
     '("descitem" 1)
     '("mbf" 1))
    (TeX-run-style-hooks
     "html"
     "epsfig"
     "latex2e"
     "art12"
     "article"
     "12pt"
     "dvips"
     "include")))

