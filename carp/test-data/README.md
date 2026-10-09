# Test fixtures

## genome_micro.fasta, custom_library.fasta, tandem_repeat_library.fasta

Inputs for the `carp` pipeline test: a 200 kB genome slice and two small, already
valid user libraries, built by `make_test_libraries.py`.

## validate_tandem_input.fasta, validate_custom_input.fasta

Inputs for the `carp_validate_repeat_library` tests. Both are the golden inputs of
upstream's own test suite, taken unchanged from
`tests/fixtures/validate_repeat_library/` of CARP 1.9.3
(https://github.com/kavonrtep/CARP, GPL-3.0), renamed with a `validate_` prefix so
they are not confused with the two valid libraries above:

| file | upstream name | what it is |
|---|---|---|
| `validate_tandem_input.fasta` | `tandem_input.fasta` | 22 *Arabidopsis* sequences from GenBank (CEN180, 5S and 45S rDNA, telomere, minisatellites), with the class written without the `Satellite/` prefix, as a real-world library would be. All 22 convert. |
| `validate_custom_input.fasta` | `custom_input.fasta` | 18 synthetic records, one per conversion path: RepeatMasker/Dfam names, REXdb/DANTE/DANTE_TIR spellings, a bare lineage name, Wicker codes, FASTA hygiene cases (description after the ID, gaps, `U`, duplicate ID), plus three records that cannot be converted (`tRNA`, no class, a protein sequence). |

Using upstream's inputs means the expected values asserted in the tool's `<tests>`
are the behaviour upstream itself pins, rather than a second interpretation of the
rules. Verified against the 1.9.3 image: the outputs this wrapper produces have the
same byte counts as upstream's expected files (35,773 B tandem FASTA, 894 B dropped
custom FASTA, 0 B for the refused library, 2,791 and 3,393 B for the two TSVs).

To refresh after an upstream change, replace `1.9.3` below:

```
B=https://raw.githubusercontent.com/kavonrtep/CARP/1.9.3/tests/fixtures/validate_repeat_library
curl -sL "$B/tandem_input.fasta" -o validate_tandem_input.fasta
curl -sL "$B/custom_input.fasta" -o validate_custom_input.fasta
```
