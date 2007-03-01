(TeX-add-style-hook "usrman"
 (lambda ()
    (LaTeX-add-bibitems
     "CoadYourdon"
     "DP"
     "ke-ri"
     "oofem"
     "pat"
     "pat2"
     "pat3"
     "pat4"
     "c++"
     "zimm")
    (LaTeX-add-labels
     "engngNummetsection"
     "Engngmodelexample"
     "materialEleemntFrame"
     "coyo"
     "genstructfig"
     "engngNummet1fig"
     "engngNummet2fig"
     "materelementFrame1"
     "microplaneFig")
    (TeX-add-symbols
     '("del" 2)
     '("mbf" 1)
     '("file" 1)
     '("attribute" 1)
     '("service" 1)
     '("class" 1)
     "refman"
     "mbff"
     "mbfx")
    (TeX-run-style-hooks
     "fancyhdr"
     "a4"
     "epsf"
     "html"
     "graphicx"
     "latex2e"
     "art12"
     "article"
     "12pt"
     "draft"
     "include")))

