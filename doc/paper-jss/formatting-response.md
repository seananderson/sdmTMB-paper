We have made all the requested changes as described below:

## Replication Material:

> The file replication-code-pcod.R contains
>
> ```
> # from remotes::install_github('seananderson/ggsidekick')
> theme_sleek <- function(base_size = 11, base_family = "") {
> ...
> }
> ```
>
> but replication-code-timings still contains
>
> ```
> ggsidekick::theme_sleek() +
> ```
>
> requiring to install the package from GitHub.

We have corrected this and no longer use `ggsidekick::`.

> The following note is obtained
>
> ```
> 2: `aes_string()` was deprecated in ggplot2 3.0.0.
> Please use tidy evaluation ideoms with `aes()`
> ```
>
> It would seem preferable to use code which avoids such notes for current
> versions of R packages.

We have corrected this to avoid this note.

## Manuscript style comments:

> Please only introduce an abbreviation within the abstract if it
> is needed again within the abstract. Otherwise, please introduce
> it within the body of the manuscript.

We have removed the GLMM abbreviation in the abstract.

> For the code layout in R publications, we typically distinguish input/output
> using Sinput/Soutput (or equivalently CodeInput/CodeOutput). Unless there are
> special reasons to format it differently, the input should use the text width
> (up to 76 or 77 characters) and be indented by two spaces, e.g.,
>
> ```
> begin{Sinput}
> R> example_model <- lm(response ~ variable1 + variable2 + variable3,
> +    weights = w, data = mydata)
> \end{Sinput}
> ```

We have corrected this and now use `R>` and `+  `.

> All captions should appear below the corresponding figure/table. The captions
> should be in sentence style and end with a period.  No additional formatting
> (such as `\emph`, f or `\it`) should be used for the caption.

We have corrected this. Fig. 1 and Table 1 captions include some `\code` and equation formatting, but we have checked with the editorial office on these.

> All table row/column headers should also be in sentence style. There should
> not be further footnote-style annotations in tables; these should all be
> placed in the caption.

We have corrected this. All footnotes are now included in the Table 1 caption.

> As a reminder, please make sure that:
> - `\proglang`, `\pkg` and `\code` have been used for highlighting throughout the
> paper (including titles and references), except where explicitly escaped.

We have gone through and made sure that all packages use `\pkg` including in section headers.

## References:

> - Journal of the Royal Statistical Society B (not: Journal of the Royal
>   Statistical Society, Series B) for Large Spatial Data Sets.” Journal of the
>   Royal Statistical Society Series B for Large Spatial Data Sets.” Journal of
>   the Royal Statistical Society Series B
> - Springer-Verlag (not: Springer) Diggle PJ, Ribeiro PJ (2007). Model-Based
>   Geostatistics. Springer. ISBN 978-0-387-48536-2.
> - John Wiley & Sons (not: Wiley, John Wiley & Sons Inc.) Wiley
> - Please make sure that all software packages are `\cite{}`'d properly.
> - All references should be in title style.

We have made all these corrections.

We kept track of all these changes in this GitHub issue:

<https://github.com/seananderson/sdmTMB-paper/issues/13>
