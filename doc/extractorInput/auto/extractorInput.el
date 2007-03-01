(TeX-add-style-hook "extractorInput"
 (lambda ()
    (TeX-add-symbols
     '("excommand" 1)
     '("mbf" 1))
    (TeX-run-style-hooks
     "latex2"
     "art10"
     "article"
     "epsf"
     "a4"
     "html"
     "include")))

