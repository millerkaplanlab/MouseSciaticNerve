#' View MouseSciaticNerve data in the scClustViz Shiny app
#'
#' A wrapper function to view one of the six \code{MouseSciaticNerve} datasets / analyses in the
#' \code{scClustViz} Shiny app.
#'
#' @param dataset One of "Inj9dBeads", "Inj9dBeadsMesenchymal", "Inj9dFACS", "InjCombined", "InjUninjMesenchymalCombined", "UninjMesenchymal".
#'   Referring to the different datasets in this package, which each represent single-cell RNAseq data collected from different Injured/Uninjured Sciatic Nerve.
#'
#' @param outPath Default = "./" (the working directory). Specify the directory
#'   used to save/load any analysis files you generate while exploring the
#'   \code{MouseSciaticNerve} data.
#'
#' @return The function causes the scClustViz Shiny GUI app to open in a
#'   seperate window.
#'
#' @examples
#'   viewMouseSciaticNerve("Inj9dBeads")
#'
#' @seealso \url{https://baderlab.github.io/scClustViz} for information on
#'   \code{scClustViz}.
#'
#' @export

viewMouseSciaticNerve <- function(dataset,outPath="./") {
  if (!dataset %in% c("Inj9dBeads", "Inj9dBeadsMesenchymal", "Inj9dFACS", "InjCombined", "InjUninjMesenchymalCombined", "UninjMesenchymal")) {
    stop('dataset must be one of "Inj9dBeads", "Inj9dBeadsMesenchymal", "Inj9dFACS", "InjCombined", "InjUninjMesenchymalCombined", "UninjMesenchymal"')
  }
  filePath <- system.file(paste0(dataset,"/",dataset,".RData"),package="MouseSciaticNerve")
  cellMarkers <- list()


  if (require("org.Mm.eg.db",quietly=T)) {
    annotationDB <- org.Mm.eg.db
    scClustViz::runShiny(filePath=filePath,
                         outPath=outPath,
                         cellMarkers=cellMarkers,
                         annotationDB=annotationDB)

  } else {
    scClustViz::runShiny(filePath=filePath,
                         outPath=outPath,
                         cellMarkers=cellMarkers)
  }
}
