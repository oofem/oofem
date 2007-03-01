(TeX-add-style-hook "oofemInput"
 (lambda ()
    (LaTeX-add-labels
     "_running_the_code"
     "_OutputFileRecord"
     "_JobDescriptionRecord"
     "_AnalysisRecord"
     "LinearStatic"
     "LinearStability"
     "EigenValueDynamic"
     "NlDEIDynamic"
     "DEIDynamic"
     "DIIDynamic"
     "IncrementalLinearStatic"
     "NonLinearStatic"
     "AdaptiveLinearStatic"
     "StationaryTransport"
     "LinearTransientTransport"
     "TransientTransport"
     "cbsIncomp"
     "staggeredproblem"
     "sparselinsolver"
     "linsolvstoragecompattable"
     "sparsesolverparams"
     "precondtable"
     "errorestimators"
     "eetypestable"
     "meshpackages"
     "ExportModulesSec"
     "_DomainRecord"
     "_OutputManagerRecord"
     "_ComponentsSizeRecord"
     "_NodeElementSideRecords"
     "_ElementsRecords"
     "_CrossSectionRecords"
     "_MaterialTypeRecords"
     "_NonlocalBarrierRecords"
     "_LoadBoundaryInitialConditions"
     "_InitialConditions"
     "_TimeFunctionsRecords"
     "ex01"
     "ex02"
     "nodecut-ex01"
     "elmentcut-ex02"
     "nodecut"
     "nodecut-lm"
     "nodecut-nlm"
     "elmentcut"
     "elmentcut-lm")
    (TeX-add-symbols
     '("PoptField" 2)
     '("Pfield" 2)
     '("PentKeywordInst" 1)
     '("PentKeywordWithVal" 2)
     '("PentKeyword" 1)
     '("Pkeyword" 2)
     '("Pkeywordnotype" 1)
     '("PfieldVal" 2)
     '("Pparam" 1)
     '("Pmode" 1)
     '("optField" 2)
     '("field" 2)
     '("entKeywordInst" 1)
     '("entKeywordWithVal" 2)
     '("entKeyword" 1)
     '("keyword" 2)
     '("keywordnotype" 1)
     '("fieldVal" 2)
     '("param" 1)
     '("mbf" 1))
    (TeX-run-style-hooks
     "html"
     "epsf"
     "graphics"
     "latex2e"
     "art10"
     "article"
     "draft"
     "include")))

