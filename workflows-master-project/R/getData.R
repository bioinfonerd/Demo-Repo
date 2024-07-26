## A collection of scripts that download relevant data off Synapse
##
## by Artem Sokolov

#' Data retrieval
#'
#' Downloads relevant data for individual workflow vignettes, usually from Synapse
#'
#' All data is downloaded into a local data/ directory
#'
#' @param vignLabel A short label specifying which vignette to get the data for. See Details.
#' @details Use the following values for \strong{vignLabel} to designate which dataset to load:
#' \describe{
#'   \item{"dge-wrangle"}{Example output from bcbio that is used by the DGE wrangling vignette}
#'   \item{"dge-pca"}{Data for the vignette that demonstrates how to apply Principal Components
#'                  Analysis to a counts matrix from a DGE experiment}
#' }
#' @export
getData <- function( vignLabel )
{
    dir.create( "data", showWarnings = FALSE )

    if( vignLabel == "dge-wrangle" )
    {
        ## Define the destination directory on a local disk
        ## Create it if it doesn't exist
        d <- "data/dge-wrangle"
        dir.create( d, showWarnings = FALSE )

        ## Download individual files to the local directory
        synapseClient::synapseLogin( rememberMe=TRUE )
        for( f in c( "syn12184051", "syn12184052", "syn12184053", "syn12184049" ) )
            synapseClient::synGet( f, downloadLocation = d )
    }
    else if( vignLabel == "dge-pca" )
    {
        ## Define the destination directory on a local disk
        ## Create it if it doesn't exist
        d <- "data/dge-pca"
        dir.create( d, showWarnings = FALSE )

        ## Download individual files to the local directory
        synapseClient::synapseLogin( rememberMe=TRUE )
        synapseClient::synGet( "syn9952382", downloadLocation = d )
        synapseClient::synGet( "syn9952383", downloadLocation = d )
    }
    else
        stop( "Unrecognized vignette label: ", vignLabel )

    invisible(NULL)
}
