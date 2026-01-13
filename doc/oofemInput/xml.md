# XML input format

The XML input format can be used everywhere where the text-based format is used. Currently, it requires oofem to be compiled with `-DUSE_XML=1`, which will links against system-installed [pugixml](https://pugixml.org) library for parsing XML. XML input will be activated for input files having the `.xml` extension automatically.

- records (lines in the text input) are expressed as entities, fields are entity attributes;
- entities are nested;
- text between elements (`PCDATA` as it is called in XML) is only used for `<Output>...</Output>` and `<Description>...</Description>` (not anywhere else);
- there is an [experimental tool](https://github.com/eudoxos/oofem-xml/) to convert from text to XML which currently passes (almost) all unit tests; that repository also contains all tests converted to XML;
- tags and attributes are case-insensitive, but message (warning) will be shown for non-matching case;
- special characters in field names are replaced by `_` to conform with XML attribute names (thus `f(t)` becomes `f_t_`, `a/c` becomes `a_c`)
- the `_` attribute denotes object to be created, e.g. `<Analysis _="StaticStructural">` (keyword argument in the text format)

In comparison to the text format:

- XML is machine readable and writeable, using usual XML tools;
- versioned (`<oofem version="...">`), with automated upgrade paths (not yet implemented);
- number of enumerated sub-itmes (such as `nelem` in a domain) are inferred (the {ref}`ComponentsSizeRecord` is therefore completely missing);
- ordering of heterogeneous records (such as whether `<Sets>` are specified before `<Elements>`) is arbitrary;
- `id`'s of enumerated elements are optional (but checked);
- XML comments `<!-- .... -->` are supported;
- XML fragments can be included via `<xi:include href="other.xml">`;
- `#%BEGIN_CHECK% ... #%END_CHECK%` (`ErrorCheckingExportModule`) is specified as XML, using `<ErrorCheck>` group, instead of special syntax.


## Field formats

Since attributes in XML are string-only, here we describe how other types are represented:

:::{list-table}
:header-rows: 1

- - type
  - description
  - notes
- - integers, floats
  - usual ASCII representation
  -
- - flags
  - attributes without value are written as `attribute=""` in XML
  - non-empty value will produce warning about not having been processed
- - boolean
  - true (`1`, `y`, `Y` and a few others), false (`0`, `n`, `N` and others)
  - unrecognized value produces error
- - 1D arrays (ints, floats)
  - space-separated
  -
- - 2D matrices
  - colon-separated rows, space-separated elements
  - e.g. `a11 a21; a12 a22`
- - dictionaries
  - colon-separated items, space-separated key and value
  - e.g `key1 value; key2 value`
- - ranges
  - comma-separated items or contiguous ranges with dash (minus)
  - e.g. `1-4,7,9-11`; as alternative: `1..4,7,9..1`
- - enums
  - specified as number of keyword
  -
:::


## Input example

* `spring01.xml`:

  ```{literalinclude} ../../tests/sm/spring01.xml
  :language: xml
  ```

* `spring01-exportmodules.xml` (included via `xi:include`):
  
  ```{literalinclude} ../../tests/sm/spring01-exportmodules.xml
  :language: xml
  ```
* `spring01-errorcheck.xml` (included via `xi:include`):
  
  ```{literalinclude} ../../tests/sm/spring01-errorcheck.xml
  :language: xml
  ```

