# segmentation functions 

require(fastseg)
require(mclust)

#' segment the methylation profile
#' 
#' The function uses fastseg() function to segment the profiles and
#' gaussian mixture modelling to cluster the segments into k components. Each
#' component ideally indicates quantitative classification of segments. Such
#' as high or low methylated regions.
#' 
#' @param score.gr \code{\link[GenomicRanges]{GRanges}} objects containing methylation scores in meta-columns
#' @param diagnostic.plot if TRUE a diagnostic plot is plotted. The plot shows
#'        metthylation and length statistics per segment group. In addition, it 
#'        shows diagnostics from mixture modeling: the density function estimated 
#'        and BIC criterion used to decide the optimum number of components
#'        in mixture modeling.
#' @param obj \code{\link[GenomicRanges]{GRanges}}, \code{\link[methylKit]{methylRaw}} or \code{\link[methylKit]{methylDiff}} object to be segmented
#' @param ... arguments to fastseg() function in fastseg package, or to densityFind
#'        in Mclust package, could be used to fine tune the segmentation algorithm.
#'        E.g. Increasing "alpha" will give more segments. 
#'        Increasing "cyberWeight" will give also more segments.
#'        For more details see fastseg and Mclust documentation.    
#'        
#' @usage methSeg(pathToFile, diagnostic.plot=TRUE, ...)
#'               
#' @return A \code{\link[GenomicRanges]{GRanges}} object with segment classification and information. 
#'        'seg.mean' column shows the mean methylation per segment.
#'        'seg.group' column shows the segment groups obtained by mixture modeling
#'        
#'        
#' @section Details:      
#'        To be sure that the algorithm will work on your data, the object should have at least 5.000 rows
#' 
#' @export
#' @docType methods
#' @rdname methSeg-methods        
methSeg<-function(obj, diagnostic.plot=TRUE, ...){
  require(fastseg)
  require(mclust)
  
  dots <- list(...)  
  
  if(class(obj)=="methylRaw"){
    obj= as(obj,"GRanges")
    mcols(obj)=100*obj$numCs/obj$coverage
  }else if(class(obj)=="methylDiff"){
    obj = as(obj,"GRanges")
    obj = sort(obj[,-1])
  }else if (class(obj) != "GRanges"){
    stop("only methylRaw, methylDiff or GRanges objects can be used in this function")
  }
  
  # match argument names to fastseg arguments
  args.fastseg=dots[names(dots) %in% names(formals(fastseg)[-1] ) ]  
  
  # match argument names to Mclust
  args.Mclust=dots[names(dots) %in% names(formals(Mclust)[-1])  ]
  
  args.fastseg[["x"]]=obj
  
  # do the segmentation
  seg.res <- do.call("fastseg", args.fastseg)
  
  # decide on number of components/groups
  args.Mclust[["score.gr"]]=seg.res
  args.Mclust[["diagnostic.plot"]]=diagnostic.plot
  dens=do.call("densityFind", args.Mclust  )
  
  # add components/group ids 
  mcols(seg.res)$seg.group=as.character(dens$classification)
  
  seg.res
}

# not needed
.methSeg<-function(score.gr,score.cols=NULL,...){
  require(fastseg)
  
  
  if(!is.null(score.cols)){
    values(score.gr)=score.gr[,score.cols]
  }
  
  seg.res <- fastseg(score.gr,...)
  
}

# finds segment groups using mixture modeling
densityFind<-function(score.gr,diagnostic.plot=T,...){
  dens = densityMclust(score.gr$seg.mean,... )
  
  if(diagnostic.plot){
    diagPlot(dens,score.gr)
  }
  dens
}


# diagnostic plot, useful for parameter trials
diagPlot<-function(dens,score.gr){
  
  scores=score.gr$seg.mean
  par(mfrow=c(2,3))
  boxplot(
    lapply(1:dens$G,function(x) scores[dens$classification==x] ),
    horizontal=T,main="methylation per group",xlab="methylation")
  
  boxplot(
    lapply(1:dens$G,function(x) log10(width(score.gr)[dens$classification==x]) ),
    horizontal=T,main="segment length per group",
    xlab="log10(length) in bp ",outline=FALSE)
  
  
  #lapply(1:dens$G,function(x) mean(width(score.gr)[dens$classification==x] ))
  
  barplot(table(dens$classification),xlab="segment groups",
          ylab="number of segments")
  plot(dens,what="density")  
  plot(dens,what="BIC")  
  
  par(mfrow=c(1,1))
  
}

#' Export segments to BED files
#' 
#' @param segments \code{\link[GenomicRanges]{GRanges}} object with segment classification and information
#' @param trackLine UCSC browser trackline
#' @param filename name of the output data
#' 
#' @return A BED files with the segmented data
#' which can be visualized in the UCSC browser 
#'
#'    
#' @export
#' @docType methods
#' @rdname methSeg2bed-methods 
methSeg2bed<-function(segments,
                      trackLine="track name='meth segments' description='meth segments' itemRgb=On",
                      filename="data/H1.chr21.chr22.trial.seg.bed"){
  require(rtracklayer)
  ramp <- colorRamp(c("gray","green", "darkgreen"))
  mcols(segments)$name=as.character(segments$seg.group)
  
  #sscores=(segments$seg.mean-min(segments$seg.mean))/(max(segments$seg.mean))-(min(segments$seg.mean))
  scores=(segments$seg.mean-min(segments$seg.mean))/(max(segments$seg.mean)-min(segments$seg.mean))
  
  mcols(segments)$itemRgb= rgb(ramp(segments$seg.mean), max = 255) 
  strand(segments)="."
  score(segments)=segments$seg.mean*100 
  
  if(is.null(trackline)){
    
    export.bed(res,filename)
  }else{
    export.bed(res,filename,
               trackLine=as(trackLine, "BasicTrackLine"))
  }
}
