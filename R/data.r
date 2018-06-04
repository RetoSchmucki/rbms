#' Values of the European portion of Metzger's high-resolution bioclimate map of the world.
#'
#' Aim: To develop a novel global spatial framework for the integration and analysis of ecological and environmental data. 
#' Location: The global land surface excluding Antarctica. Methods: A broad set of climate-related variables were considered 
#' for inclusion in a quantitative model, which partitions geographic space into bioclimate regions. Statistical screening 
#' produced a subset of relevant bioclimate variables, which were further compacted into fewer independent dimensions using 
#' principal components analysis (PCA). An ISODATA clustering routine was then used to classify the principal components into 
#' relatively homogeneous environmental strata. The strata were aggregated into global environmental zones based on the attribute 
#' distances between strata to provide structure and support a consistent nomenclature. Results: The global environmental 
#' stratification (GEnS) consists of 125 strata, which have been aggregated into 18 global environmental zones. The stratification 
#' has a 30 arcsec resolution (equivalent to 0.86 km2 at the equator). Aggregations of the strata were compared with nine existing 
#' global, continental and national bioclimate and ecosystem classifications using the Kappa statistic. Values range between 
#' 0.54 and 0.72, indicating good agreement in bioclimate and ecosystem patterns between existing maps and the GEnS. Main conclusions: 
#' The GEnS provides a robust spatial analytical framework for the aggregation of local observations, identification of gaps in 
#' current monitoring efforts and systematic design of complementary and new monitoring and research.
#' The original dataset is available for non-commercial use through the GEO portal (http://www.geoportal.org). 
#' Value's classification key is available in \code{metzger_v3_class} dataset.
#'
#' @format values used to build a raster object \code{raster_data/metzger_v3_europe.tif} using the empty raster \code{metger_v3_europe.rda}:
#' \describe{
#'   \item{value}{Integer associated classification is found in the \code{metzger_v3_class} dataset}
#' }
#' @references Metzger, M.J., Bunce, R.G.H., Jongman, R.H.G., Sayre, R., Trabucco, A. & Zomer, R. (2013) A high-resolution bioclimate map of the world: a unifying framework 
#' for global biodiversity research and monitoring. Global Ecology and Biogeography, 22, 630–638.
#' @source Global map is available from \code{\url{https://edinburgh-innovations.ed.ac.uk/project/bioclimate-world-map}}.

"metzger_v3_europe_values"