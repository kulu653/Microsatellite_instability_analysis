score_easy_kul_A27_125bp_to_150_with_shutter<-
  function (my.inds, cols = 1, n.inds = NULL, panel = NULL, shift = 0.8,ladder, channel.ladder = NULL, ploidy = 2, left.cond = c(0.6, 3), right.cond = 0.35, warn = FALSE, window = 0.5, init.thresh = 200,ladd.init.thresh = 200, method = "iter2", env = parent.frame(), 
            plotting = TRUE, electro = FALSE, pref = 4)
  {
    
    oldw <- getOption("warn")
    options(warn = -1)
    ci.upp=1.96; ci.low=1.96; dev=50 ;thresh=NULL
    if(length(n.inds) > length(my.inds)){
      print(paste("Hey! you are trying to examine more individuals than the ones you actually read? You selected in 'n.inds' argument", length(n.inds), "individuals but you only provided", length(my.inds), " individuals. Please select a number of individuals smaller or same size than the ones contained in 'my.inds' argument"))
      stop
    }else{
      cat(paste("\n1) You have used a shift of", shift, "base pairs. All peaks at that distance from the tallest peak will be ignored and be considered noise. \n2) In addition the window used is", window, ". Which means that all peaks closer by that distance to panel peaks will be accounted as peaks. \n3) Remember using the get.scores() function to extract the results from this output. \n\n"))
      #print(paste("In addition the window used is", window, ". Which means that all peaks closer by that distance to panel peaks will be accounted as peaks"))
    }
    if(method == "ci"){
      print(paste("Please make sure you have used the same 'dev' value you found convenient for your ladder detection or probably your call will not match"))
    }
    # the shift is used to declare when peaks will be considered the same peak and only the tallest will be picked
    # window willl be used to declase a peak the same for each base pair, i.e 170.4 is 170 same than 169.7
    ### makes a subset with the desired plants
    
    ## channel where the ladder is located
    if(is.null(channel.ladder )){
      channel.ladder <- dim(my.inds[[1]])[2]
    }else{ channel.ladder <- channel.ladder}
    
    if(dim(my.inds[[1]])[2] < channel.ladder){
      print(paste("ERROR MY FRIEND!! you have indicated an argument channel.ladder=5, but your data contains less channels/colors"))
      stop
    }
    ## number of samples to do
    if(is.null(n.inds )){
      n.inds <- c(1:length(my.inds))
    }else{n.inds <- n.inds}
    ## thresh
    if(is.null(thresh)){
      thresh <- rep(list(c(1,1,1,1,1)), length(my.inds))
    }else{thresh <- thresh}
    ####################
    ## initialize the progress bar
    count <- 0
    tot <- length(n.inds)
    pb <- txtProgressBar(style = 3)
    setTxtProgressBar(pb, 0)
    #####################
    
    my.inds2 <- list(NA)
    thresh2 <- list(NA)
    for (i in 1:length(n.inds)) {
      ###################
      count <- count + 1
      ##################
      v1 <- n.inds[i]
      my.inds2[[i]] <- my.inds[[v1]]
      names(my.inds2)[i] <- names(my.inds)[i]
      #################################
      setTxtProgressBar(pb, (count/tot)*.25)### keep filling the progress bar
      ################################
    }
    #my.inds <- my.inds2
    #if (!require("zoom")) {
    # install.packages("zoom")
    #  require("zoom")
    #}
    ncfp <- c("COLOR 1", "COLOR 2", "COLOR 3", "COLOR 4", "COLOR 5", "COLOR 6")
    cfp <- c("cornflowerblue", "chartreuse4", "gold2", "red", "orange", "purple")
    col.list <- list(NA)
    att1 <- numeric()
    ######################################################################################################
    ### this part extracts all the models for each single plant in my plants
    list.data <- list(NA)
    if(exists("list.data.covarrubias")){
      list.data <- env$list.data.covarrubias
    }else{
      list.ladders <- lapply(my.inds2, function(x){y <- x[,channel.ladder]; return(y)})
      # extract ladder channels for all plants
      list.data <- lapply(list.ladders, find.ladder, ladder=ladder, ci.upp=ci.upp, ci.low=ci.low, draw=F, dev=dev, warn=warn, method=method,init.thresh=ladd.init.thresh)
    }
    # this models uses indexes and predicts base airs
    list.models <- lapply(list.data, function(da){y <- da[[3]]; x <- da[[1]];mod <- lm(y~ I(x) + I(x^2) + I(x^3) + I(x^4) + I(x^5), data=da); return(mod)})
    # this models uses pairs and predicts indexes
    list.models.inv <- lapply(list.data, function(da){x <- da[[3]]; y <- da[[1]];mod <- lm(y~ x, data=da); return(mod)})
    #########################################################################
    xx <- lapply(my.inds2, function(x, cols){1:length(x[,cols])}, cols=cols)
    newxx <- numeric()
    newyy <- numeric()
    new.whole.data <- list(NA)
    for(h in 1:length(xx)){
      h1 <- n.inds[h]
      ###################
      count <- count + 1
      ###################
      # h1 is used to make sure it matches if nplants is a weird vector
      newxx <- as.vector(predict(list.models[[h1]], newdata=data.frame(x=xx[[h]])))
      newyy <- my.inds2[[h]][,cols]
      new.whole.data[[h]] <- list(xx=newxx, yy=newyy)
      ################
      setTxtProgressBar(pb, (count/tot)*.25)### keep filling the progress bar
      ################################
    }
    
    top <- max(unlist(lapply(new.whole.data, function(x){max(x$yy)})))
    bott <- min(unlist(lapply(new.whole.data, function(x){min(x$yy)})))
    ### STARTS REDUCTION OF THE DATA ###
    list.weis <- list(NA) 
    lower.bounds <- numeric()
    for(k in 1:length(my.inds2)){
      #sh <- as.vector(predict(list.models.inv[[k]], newdata=data.frame(x=shift))) - list.models.inv[[k]]$coefficients[[1]]
      #if(length(panel) > 0){
      #pani <- c(min(panel)-1,  max(panel)+1)
      #maxpeak.in.pani <- max(new.whole.data[[k]][[2]][which(new.whole.data[[k]][[1]] > pani[1] & new.whole.data[[k]][[1]] < pani[2])])
      ## just create newthrsholds for plants who actually have peaks in that panel region
      #if( maxpeak.in.pani > init.thresh){
      # newtt <- threshs(my.plant=new.whole.data[[k]], min.thre=init.thresh, panel=pani, ci=thresh[[k]][cols])  
      #}else{newtt <- init.thresh}
      newtt <- init.thresh
      #####################################
      lower.bounds[k] <- newtt
      #}else{newtt <- init.thresh;  lower.bounds[k] <- newtt}
      plant <- big.peaks.col(new.whole.data[[k]]$yy, newtt)
      plant$wei <- new.whole.data[[k]]$xx[plant$pos]
      plant <- separate(plant, shift, type="bp")
      
      #plant$wei <- new.whole.data[[k]]$xx[plant$pos]
      list.weis[[k]] <- plant
      
      
      ###################
      if(plotting == TRUE){
        count <- count + 1
      }else{count <- count + 2}
      
      setTxtProgressBar(pb, (count/tot)*.25)### keep filling the progress bar
      
      ################################
    }
    # round digits
    list.weis <- lapply(list.weis, function(x){x$wei <- round(x$wei,digits=4); return(x)})
    names(list.weis) <- names(my.inds2)
    
    if(length(panel) > 0){
      list.weis <- lapply(list.weis, reals, panel=panel, shi=shift, ploidy=ploidy, left.cond=left.cond, right.cond=right.cond, window=window)
      list.weis2 <- lapply(list.weis, FUN=homo.panel, panel=panel, window=window)
      
      
      
    }else{
      list.weis2 <- list.weis
      
    }
    
    
    # END OF THE ANALYSIS
    
    Plot_name<- basename(getwd())
    if(plotting == TRUE){
      #### PLOTTING PART ###
      
      #png(file='plot2.png')
      pdf(file = file.path("C:/Users/kshresth/MSI_analysis/Results/electrograms with shutters_pdfs",paste (Plot_name,cfp[cols],".pdf")), paper="a4",height= 12, width =8)
      par(mfrow=c(4,2))
      if(length(panel) > 0){
        
        xm <- round(min(panel, na.rm = TRUE)-10, digits=0)
        xl <- round(max(panel, na.rm = TRUE)+10, digits=0)
      }else{xm <- 0; xl <- max(ladder)}
      for(g in 1:length(n.inds)){
        hh4 <- n.inds[g]
        
        if(length(which(new.whole.data[[g]]$xx > xm & new.whole.data[[g]]$xx < xl)) > 0){
          mylim <- max(new.whole.data[[g]]$yy[which(new.whole.data[[g]]$xx > xm & new.whole.data[[g]]$xx < xl)], na.rm = TRUE) +100
        }else{mylim=1000}
        
        if(is.infinite(mylim)){mylim=1000}
        plot( new.whole.data[[g]]$xx, new.whole.data[[g]]$yy, type="l", col=cfp[cols], xaxt="n",
              xlim=c(125,150),  ylim=c(0,4000), ylab="Intensity",  # real code ma 0 ko sato ylim=c(lower.bounds[g],mylim +100)
              xlab=names(list.models)[hh4], lwd=2, las=2)
        axis(1,at=c(xm:xl), labels=xm:xl, cex.axis=0.8)
        #rect(xleft=(list.weis2[[g]]$wei-window),ybottom=(bott-200), xright=(list.weis2[[g]]$wei+window), ytop=(top+1000), col=transp("lightpink",0.3), border=NA)
        abline(v=list.weis[[g]]$wei, lty=3, col="blue", cex=0.5)
        #MERO EDIR FOR SIZE SHIFT
        text(x=list.weis[[g]]$wei- 0.5, y=600, srt= 90, labels=round(list.weis[[g]]$wei, digits = 2),cex=0.75)
        ###################
        
        #########
        count <- count + 1
        setTxtProgressBar(pb, (count/tot)*.25)### keep filling the progress bar
        
       
      }
      
      close(pb) # close the progress bar
      options(warn = oldw)
      return(list.weis2)
    }
    
  }

