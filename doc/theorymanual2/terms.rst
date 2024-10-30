Term documentation (MPM module)
###############################


+-----------------+------------------------------------------------------------------+------------------------------------------------------------------+
| Name [class]    | Definition                                                       | Parameters/Notes                                                 |
+=================+==================================================================+==================================================================+
| ``BTSigTerm``   | :math:`\int_\Omega \nabla^s w\ \sigma(\nabla^s u)`               | :math:`w` - test variable, :math:`u` -  unknown variable,        |
|                 |                                                                  | The :math:`\sigma` is (nonlinear) operator evaluated by a        |
|                 |                                                                  | constitutive model.                                              |
|                 |                                                                  | Supported constitutive modes: ``_3dMat, _3dUP, _2dUP``.          |
+-----------------+------------------------------------------------------------------+------------------------------------------------------------------+
| ``gNTfTerm``    | :math:`\int_\Omega \nabla w\ \mathbf{f}(\nabla p)`               | :math:`w` - test variable, :math:`p` -  unknown variable,        |
|                 |                                                                  | The :math:`\bf{f}` is (nonlinear) constitutive operator          |
|                 |                                                                  | evaluated by a constitutive model.                               |
|                 |                                                                  |                                                                  |
+-----------------+------------------------------------------------------------------+------------------------------------------------------------------+
| ``BTamNTerm``   | :math:`\int_\Omega \nabla^s w\ \alpha\mathbf{m}\ p`              | :math:`w` - test variable, :math:`p` -  unknown variable,        |
|                 |                                                                  | The :math:`\alpha` is constitutive parameter provided by material|
|                 |                                                                  | model, The vector :math:`\mathbf{m}=\{1,1,1,0,0,0\}^T`           |
+-----------------+------------------------------------------------------------------+------------------------------------------------------------------+
| ``NTamTBTerm``  | :math:`\int_\Omega w\ \alpha\mathbf{m}\ \nabla^s \dot{u}`        | :math:`w` - test variable, :math:`u` -  unknown variable,        |
|                 |                                                                  | The :math:`\alpha` is constitutive parameter provided by material|
|                 |                                                                  | model, The vector :math:`\mathbf{m}=\{1,1,1,0,0,0\}^T`           |
+-----------------+------------------------------------------------------------------+------------------------------------------------------------------+
