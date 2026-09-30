Output and job description records
==================================

The first two records of every input file carry no keyword: they are recognised
by their position alone.

.. _OutputFileRecord:

Output file record
------------------

A character string giving the path to the output file.  If a file of the same
name already exists, it is overwritten.

In parallel mode the partition name is **not** appended to this path, unlike the
input file name given with ``-f`` (see :ref:`running-the-code`).  Each partition
input file should therefore name its own output file.

.. _JobDescriptionRecord:

Job description record
----------------------

A character string describing the job.  The description is reproduced in the
output file.
