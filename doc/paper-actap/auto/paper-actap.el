(TeX-add-style-hook "paper-actap"
 (lambda ()
    (LaTeX-add-bibitems
     "CoadYourdon"
     "ke-ri"
     "pat"
     "c++"
     "zimm")
    (LaTeX-add-labels
     "engngNummetsection"
     "materialEleemntFrame"
     "coyo"
     "genstructfig"
     "engngNummet1fig"
     "engngNummet2fig"
     "materelementFrame1"
     "microplaneFig")
    (TeX-add-symbols
     '("service" 1)
     '("class" 1))
    (TeX-run-style-hooks
     "latex2"
     "art12"
     "article"
     "12pt"
     "epsf"
     "html"
     "a4"
     "include")))

