# nocov start

cite = function(entry){

 out <- switch(
  entry,
  breiman_2001 = utils::bibentry(
   "article",
   title       = "Random Forests",
   author      = "Breiman, Leo",
   year        = "2001",
   journal     = "Machine Learning",
   volume      = "45",
   number      = "1",
   pages       = "5--32",
   # doi         = "10.1023/A:1010933404324",
   issn        = "1573-0565"
  ),
  ishwaran_2008 = utils::bibentry(
   "article",
   # doi         = "10.1214/08-aoas169",
   year        = "2008",
   month       = "9",
   publisher   = "Institute of Mathematical Statistics",
   volume      = "2",
   number      = "3",
   author      = "Hemant Ishwaran and Udaya B. Kogalur and Eugene H. Blackstone and Michael S. Lauer",
   title       = "Random survival forests",
   journal     = "The Annals of Applied Statistics"
  ),
  jaeger_2019 = utils::bibentry(
   "article",
   # doi           = "10.1214/19-aoas1261",
   year          = "2019",
   month         = "9",
   publisher     = "Institute of Mathematical Statistics",
   volume        = "13",
   number        = "3",
   author        = "Byron C. Jaeger and D. Leann Long and Dustin M. Long and Mario Sims and Jeff M. Szychowski and Yuan-I Min and Leslie A. Mcclure and George Howard and Noah Simon",
   title         = "Oblique random survival forests",
   journal       = "The Annals of Applied Statistics"
  ),
  jaeger_2022 = utils::bibentry(
   "article",
   title         = "Accelerated and interpretable oblique random survival forests",
   author        = "Byron C. Jaeger and Sawyer Welden and Kristin Lenoir and Jaime L. Speiser and Matthew W. Segar and Ambarish Pandey and Nicholas M. Pajewski",
   journal       = "Journal of Computational and Graphical Statistics",
   # doi           = "10.1080/10618600.2023.2231048",
   year          = "2023",
   month         = "8",
   publisher     = "Taylor & Francis",
   pages         = "1--16"
  ),
  hooker_2021 = utils::bibentry(
   "article",
   title     = "Unrestricted permutation forces extrapolation: variable importance requires at least one more model, or there is no free variable importance",
   author    = "Hooker, Giles and Mentch, Lucas and Zhou, Siyu",
   journal   = "Statistics and Computing",
   volume    = "31",
   pages     = "1--16",
   year      = "2021",
   publisher = "Springer"
  ),
  harrell_1982 = utils::bibentry(
   "article",
   title     = "Evaluating the yield of medical tests",
   author    = "Harrell, Frank E and Califf, Robert M and Pryor, David B and Lee, Kerry L and Rosati, Robert A",
   journal   = "Jama",
   volume    = "247",
   number    = "18",
   pages     = "2543--2546",
   year      = "1982",
   publisher ="American Medical Association"
  ),
  menze_2011 = utils::bibentry(
   "inproceedings",
   title        = "On oblique random forests",
   author       = "Menze, Bjoern H and Kelm, B Michael and Splitthoff, Daniel N and Koethe, Ullrich and Hamprecht, Fred A",
   booktitle    = "Machine Learning and Knowledge Discovery in Databases: European Conference, ECML PKDD 2011, Athens, Greece, September 5-9, 2011, Proceedings, Part II 22",
   pages        = "453--469",
   year         = "2011",
   organization = "Springer"
  ),
  simon_2011 = utils::bibentry(
   "article",
   title     = "Regularization paths for Cox's proportional hazards model via coordinate descent",
   author    = "Simon, Noah and Friedman, Jerome and Hastie, Trevor and Tibshirani, Rob",
   journal   = "Journal of statistical software",
   volume    = "39",
   number    = "5",
   pages     = "1",
   year      = "2011",
   publisher = "NIH Public Access"
  ),
  horst_2022 = utils::bibentry(
   "article",
   title   = "Palmer Archipelago Penguins Data in the palmerpenguins R Package-An Alternative to Anderson's Irises",
   author  = "Horst, Allison M and Hill, Alison Presmanes and Gorman, Kristen B",
   journal = "R Journal",
   volume  = "14",
   number  = "1",
   year    = "2022"
  ),
  penguins_2020 = utils::bibentry(
   "manual",
   title  = "palmerpenguins: Palmer Archipelago (Antarctica) penguin data",
   author = "Allison Marie Horst and Alison Presmanes Hill and Kristen B Gorman",
   year   = "2020",
   note   = "R package version 0.1.0",
   url    = "https://allisonhorst.github.io/palmerpenguins/"
  ),
  greenwell_2018 = utils::bibentry(
   "article",
   title   = "A simple and effective model-based variable importance measure",
   author  = "Greenwell, Brandon M and Boehmke, Bradley C and McCarthy, Andrew J",
   journal = "arXiv preprint arXiv:1805.04755",
   year    = "2018"
  ),
  stop("unrecognized entry")
 )

 format(out)

}

# nocov end
