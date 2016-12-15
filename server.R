
#
# This is a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
  # step 1 callibration:
    # input table with calibration and serial dillutions
    # output table with linear trend and index for regression model
  # step 2 quantification:
    # input model index
    # output quantification



options(shiny.maxRequestSize=1000*1024^2) 
library(shiny)
library(sendmailR)
library(chemCal) 

shinyServer(function(input, output) {
  
  
  dataInput <- reactive({
    if(length(input$file1)) {
        inFile <- input$file1
        if(!input$gnps) {
          opt2 <- read.csv(inFile$datapath)  
        } else {
          optmp <- read.delim(inFile$datapath)  
          opt2 <- t(optmp[,-1])
          colnames(opt2) <- sub("-", "\\.", optmp[,1])
          
          inFile11 <- input$file11
          conc <- read.csv(inFile11$datapath)  
          
          if(rownames(opt2)[nrow(opt2)]=="X"){ 
            opt2 <- data.frame(SampleID=rownames(opt2)[-nrow(opt2)], conc, opt2[-nrow(opt2),])
          } else {
            opt2 <- data.frame(SampleID=rownames(opt2), conc, opt2)
           }
          rownames(opt2) <- NULL
        }
        
      mzw <- input$mz1
      rtw <- input$rt1
        
     match.peaks <- function(query, ref, mz.tol = mzw, rt.tol = rtw) {
        mzdiff <- (abs(ref[,1]-query[1])/query[1])*10^6
        rtdiff <- abs(ref[,2]-query[2])
        mtch <- which((mzdiff <= mz.tol) & (rtdiff <= rt.tol))
        # gives an unique match
        if(length(mtch)>1) mtch <- mtch[which.min(mzdiff[mtch])]
        mtch
    }
    mexp <- data.frame(suffix=c("pM", "nM", "uM"), fac=c(10^-12, 10^-9, 10^-6))

    if(length(input$file2)) {
      inFile <- input$file2
      exp <- read.csv(inFile$datapath)  
      exp2 <- as.matrix(exp[, 1:2])
    }
     
     mn <- sub("_", "", opt2[,2])
     mstep1 <- as.numeric(gsub("\\D", '', mn))
     mstep2 <- mexp[sapply(gsub("\\d", '', mn), function(x) which(mexp[,1]==x)),2]
    
     vec <- order(mstep1 * mstep2)
     
     summ <- rep("", 5) 
    
       for(j in 3:ncol(opt2)) {
         d1 <- data.frame(x=as.vector(as.matrix(log10(as.numeric(opt2[,j])))), 
				              y=log10((mstep1*mstep2)), labels=mn)
         d1 <- d1[vec,]
         x <- d1$x 
         y <- d1$y 
         y <- y[x!=-Inf]
         x <- x[x!=-Inf]
         if(length(y)>1) {
           s1 <- summary((lm(y~x)))
           summ <- rbind(summ, c(j, colnames(opt2)[j] ,  round(s1$adj.r.squared, 3), s1$coefficients[2,4], s1$coefficients[,1][2]) )
         } else{
           summ <- rbind(summ, c(j, colnames(opt2)[j],  0, 0, 0))
         }
       }
                 
     # Filter the results by R2
     summ <- summ[-1,][!is.nan(as.numeric(summ[-1,3])),]
     selsumm <- summ[as.numeric(as.matrix(summ[,3]))>0,]
     
     # Have to import header properly, to avoid R weird replacement
     name2fet <- function(fet) {
         strf <- strsplit(fet, "\\.") 
         mz <- paste(strf[[1]][1:2], collapse = ".")
         mz <- as.numeric(sub("^X", "", mz))
         rt <- as.numeric(paste(strf[[1]][3], collapse = "."))
         c(mz, rt)
      }
 
                                                          
    if(length(input$file2)) {
      
       fetm <- as.vector(as.matrix(selsumm[,2]))
       fetm <- do.call(rbind, lapply(fetm, name2fet))
      
       if(input$omass) {
           # it is only mass already, change later
           # have to have different format for exp2 here
           m.list.mz.rt_opt <- apply(exp2, 1, match.peaks, fetm, rt.tol = 1000)
       } else {
           m.list.mz.rt_opt <- apply(exp2, 1, match.peaks, fetm)
       }
       if(length(m.list.mz.rt_opt)) {
         selsumm <- cbind(selsumm, "") 
         if(input$omass) {
          for(i in 1:length(m.list.mz.rt_opt)) if(length(m.list.mz.rt_opt[[i]])) selsumm[m.list.mz.rt_opt[[i]], ncol(selsumm)] <- as.character(exp[i,2])
         } else {
          for(i in 1:length(m.list.mz.rt_opt)) if(length(m.list.mz.rt_opt[[i]])) selsumm[m.list.mz.rt_opt[[i]], ncol(selsumm)] <- as.character(exp[i,3])
         }
         selsumm[,2] <- apply(fetm, 1, function(x) paste(x, collapse = "_"))
       } else {
         selsumm <- cbind(selsumm, "") 
         selsumm[,2] <- apply(fetm, 1, function(x) paste(x, collapse = "_"))
       }
    } else {
      
       fetm <- as.vector(as.matrix(selsumm[,2]))
       fetm <- do.call(rbind, lapply(fetm, name2fet))
       
       selsumm <- cbind(selsumm, "") 
       selsumm[,2] <- apply(fetm, 1, function(x) paste(x, collapse = "_"))
      
    }
    list(selsumm=selsumm, mstep1=mstep1, mstep2=mstep2, vec=vec, opt2=opt2, mn=mn)
  }
  })
  
 
  output$data_file <- downloadHandler(
    filename = 'curve_output.zip',
    content = function(fname) {
      tmpdir <- tempdir()
      drname <- paste(sample(letters, 10), collapse="") 
      dir.create(paste0(tmpdir, "/", drname))
      setwd(paste0(tmpdir, "/", drname))
      
      dir.create("curve_output") 

      selsumm <- dataInput()$selsumm
      opt2 <- dataInput()$opt2
      vec <- dataInput()$vec
      mn <- dataInput()$mn
      mstep1 <- dataInput()$mstep1
      mstep2 <- dataInput()$mstep2
       
       pdf("curve_output/calibration_curve.pdf")
       summ <- rep("", 7)
       opar <- par()
       modelList <- list()

       for(i in 1:nrow(selsumm)) {
         # This line controls whether to plot lines that are flat 
         # as we should expect inclined lines
         # filter by R2
         if(as.numeric(selsumm[i,5]) < input$slope | as.numeric(selsumm[i,3]) < input$r2) next
         j <- as.numeric(selsumm[i,1])
      
         d1 <- data.frame(x=as.vector(as.matrix(log10(as.numeric(opt2[,j])))), 
				    y=log10((mstep1*mstep2)), labels=mn)
         d1 <- d1[vec,]

      	 x <- d1$x 
      	 y <- d1$y 
      	 y <- y[x!=-Inf]
      	 x <- x[x!=-Inf]
      	 l0 <- lm(y~x)
      
      	 modelList[[selsumm[i,1]]] <- l0 
      
      	if(input$ptype=="calplot") {
      		 calplot(l0, xlim=c(min(x), max(x)), xlab="Log intensity", ylab="Log Concentration [M]") 
      	} else {
      		 lb <- table(d1$labels)
      		 # assumes same number of technical replicates
      		 x <- rep(1:length(lb), each=lb[1])
      		 
      		 par(las=2)
      	      
      		 ym <- tapply(d1$x, d1$labels, mean)[as.character(unique(d1$labels))]
      		 plot(1:length(lb), ym, pch=19,  ylab="Log intensity", xlab="", xaxt="n") 
      		 sd <- tapply(d1$x, d1$labels, sd)[as.character(unique(d1$labels))]
      		 segments(1:length(lb), ym-sd, 1:length(lb), ym+sd)
      		 epsilon = 0.02
      		 segments(1:length(lb)-epsilon,ym-sd,1:length(lb)+epsilon,ym-sd)
      		 segments(1:length(lb)-epsilon,ym+sd,1:length(lb)+epsilon,ym+sd)
      	      
      		 axis(1, at=1:length(lb), labels=as.character(unique(d1$labels)), cex=0.8)
      	      
      		 # allow just concentrations with 3 replicates to model
      		 for(k in 1:nrow(d1)) if( d1[k,1]==-Inf) d1[d1[,3]==d1[k,3],2] <- -Inf
      	      
      		 y <- d1$x
      	      
      		 x <- x[y!=-Inf]
      		 y <- y[y!=-Inf]
      	      
      		 l1 <- lm(y~x) 
           abline(l1)
      	}

        if(selsumm[i,6]!="") {
          title(paste(selsumm[i,2], "\n",  selsumm[i,6], "\n",  
          "R2", round(as.numeric(selsumm[i,3]), 3), "p-value", format(as.numeric(selsumm[i,4]), scientific=TRUE, digits=3) ))
         } else {
          title(paste(selsumm[i,2], "\n",  
          "R2", round(as.numeric(selsumm[i,3]), 3), "p-value", format(as.numeric(selsumm[i,4]), scientific=TRUE, digits=3) ))
         }

         summ <- rbind(summ, c(selsumm[i,c(1,2,6)], round(as.numeric(selsumm[i,3]), 3), 
			                  format(as.numeric(selsumm[i,4]), scientific=TRUE, digits=3),
                               10^min(d1$x), 10^max(d1$x))
		              )
       }
       dev.off()
       cnames <- c("Feat_pos", "Mz_rt", "Name", "R2", "p-value", "MinInt", "MaxInt")

       if(!is.null(nrow(summ))) {
         summ <- matrix(summ[-1,], ncol=7)
         colnames(summ) <- cnames
         write.csv(summ, "curve_output/calibration_curve_table.csv", row.names = FALSE)
         save(modelList, file="curve_output/modelList.rda")
       } else{
         write.csv("No increasing linear trend found.", "curve_output/calibration_curve_table.csv", row.names = FALSE)
       }
      
      zip(zipfile = fname, files="curve_output", flags = "-r")
      
      
      if(file.exists(paste0(fname, ".zip"))) {file.rename(paste0(fname, ".zip"), fname)}
    },
    contentType = "application/zip"
  )   
  
  EmailOutput <- reactive({
      input$goMail 
    
      tmpdir <- tempdir()
      drname <- paste(sample(letters, 10), collapse="") 
      dir.create(paste0(tmpdir, "/", drname))
      setwd(paste0(tmpdir, "/", drname))
      
      dir.create("curve_output") 

      selsumm <- dataInput()$selsumm
      opt2 <- dataInput()$opt2
      vec <- dataInput()$vec
      mn <- dataInput()$mn
      mstep1 <- dataInput()$mstep1
      mstep2 <- dataInput()$mstep2
       
       pdf("curve_output/calibration_curve.pdf")
       summ <- rep("", 7)
       opar <- par()
       modelList <- list()

       for(i in 1:nrow(selsumm)) {
         # This line controls whether to plot lines that are flat 
         # as we should expect inclined lines
         # filter by R2
         if(as.numeric(selsumm[i,5]) < input$slope | as.numeric(selsumm[i,3]) < input$r2) next
         j <- as.numeric(selsumm[i,1])
      
         d1 <- data.frame(x=as.vector(as.matrix(log10(as.numeric(opt2[,j])))), 
				    y=log10((mstep1*mstep2)), labels=mn)
         d1 <- d1[vec,]

      	 x <- d1$x 
      	 y <- d1$y 
      	 y <- y[x!=-Inf]
      	 x <- x[x!=-Inf]
      	 l0 <- lm(y~x)
      
      	 modelList[[selsumm[i,1]]] <- l0 
      
      	if(input$ptype=="calplot") {
      		 calplot(l0, xlim=c(min(x), max(x)), xlab="Log intensity", ylab="Log Concentration [M]") 
      	} else {
      		 lb <- table(d1$labels)
      		 # assumes same number of technical replicates
      		 x <- rep(1:length(lb), each=lb[1])
      		 
      		 par(las=2)
      	      
      		 ym <- tapply(d1$x, d1$labels, mean)[as.character(unique(d1$labels))]
      		 plot(1:length(lb), ym, pch=19,  ylab="Log intensity", xlab="", xaxt="n") 
      		 sd <- tapply(d1$x, d1$labels, sd)[as.character(unique(d1$labels))]
      		 segments(1:length(lb), ym-sd, 1:length(lb), ym+sd)
      		 epsilon = 0.02
      		 segments(1:length(lb)-epsilon,ym-sd,1:length(lb)+epsilon,ym-sd)
      		 segments(1:length(lb)-epsilon,ym+sd,1:length(lb)+epsilon,ym+sd)
      	      
      		 axis(1, at=1:length(lb), labels=as.character(unique(d1$labels)), cex=0.8)
      	      
      		 # allow just concentrations with 3 replicates to model
      		 for(k in 1:nrow(d1)) if( d1[k,1]==-Inf) d1[d1[,3]==d1[k,3],2] <- -Inf
      	      
      		 y <- d1$x
      	      
      		 x <- x[y!=-Inf]
      		 y <- y[y!=-Inf]
      	      
      		 l1 <- lm(y~x) 
           abline(l1)
      	}

        if(selsumm[i,6]!="") {
          title(paste(selsumm[i,2], "\n",  selsumm[i,6], "\n",  
          "R2", round(as.numeric(selsumm[i,3]), 3), "p-value", format(as.numeric(selsumm[i,4]), scientific=TRUE, digits=3) ))
         } else {
          title(paste(selsumm[i,2], "\n",  
          "R2", round(as.numeric(selsumm[i,3]), 3), "p-value", format(as.numeric(selsumm[i,4]), scientific=TRUE, digits=3) ))
         }

         summ <- rbind(summ, c(selsumm[i,c(1,2,6)], round(as.numeric(selsumm[i,3]), 3), 
			                  format(as.numeric(selsumm[i,4]), scientific=TRUE, digits=3),
                               10^min(d1$x), 10^max(d1$x))
		              )
       }
       dev.off()
       cnames <- c("Feat_pos", "Mz_rt", "Name", "R2", "p-value", "MinInt", "MaxInt")

       if(!is.null(nrow(summ))) {
         summ <- matrix(summ[-1,], ncol=7)
         colnames(summ) <- cnames
         write.csv(summ, "curve_output/calibration_curve_table.csv", row.names = FALSE)
         save(modelList, file="curve_output/modelList.rda")
       } else{
         write.csv("No increasing linear trend found.", "curve_output/calibration_curve_table.csv", row.names = FALSE)
       }
      
      fname = 'curve_output.zip'
      zip(zipfile = fname, files="curve_output", flags = "-r")
      
      from = "<mysterious-cloud@dal.ca>"
      to = paste0("<", input$mail, ">")
      
      control=list(smtpServer="ASPMX.L.GOOGLE.COM")
      subject <- "Calibration curves"
      body <- "Here are your calibration curves, please let us know if the formatting is correct."
      #sendmail(from=from,to=to,subject=subject,msg=body,control=control)
      attachmentPath <- "curve_output.zip"
      attachmentName <- "curve_output.zip"
      attachmentObject <- mime_part(x=attachmentPath,name=attachmentName)
      bodyWithAttachment <- list(body,attachmentObject)
      
      sendmail(from=from,to=to,subject=subject,msg=bodyWithAttachment,control=control)
  })
  
  observeEvent(input$goMail, {
    EmailOutput()
  })
  
  output$data_file2 <- downloadHandler(
    filename = function() { paste("quantification", '.csv', sep='') },
    content = function(file) {
      #print("I'm here")
      
      name2fet <- function(fet) {
         strf <- strsplit(fet, "\\.") 
         mz <- paste(strf[[1]][1:2], collapse = ".")
         mz <- as.numeric(sub("^X", "", mz))
         rt <- as.numeric(paste(strf[[1]][3], collapse = "."))
         c(mz, rt)
      }
 
     
    if(length(input$file3)) {
      inFile <- input$file3
      #curve <- read.csv(inFile$datapath)  
      curve <- read.csv(unz(inFile$datapath, "curve_output/calibration_curve_table.csv"))
      load(unz(inFile$datapath, "curve_output/modelList.rda"))
      
      pcurve <- t(sapply(curve[,2], function(x) strsplit(as.character(x), "_")[[1]]))
      pcurve <- apply(pcurve, 2,  as.numeric)
      
      inFile <- input$file4
      qtn <- read.csv(inFile$datapath)  
      
       pqtn <- as.vector(colnames(qtn)[-1])
       pqtn <- do.call(rbind, lapply(pqtn, name2fet))
      
      mzw <- input$mz
      rtw <- input$rt
      
      match.peaks <- function(query, ref, mz.tol = mzw, rt.tol = rtw) {
        #mzdiff <- abs(ref[,1]-query[1])
        mzdiff <- (abs(ref[,1]-query[1])/query[1])*10^6
        rtdiff <- abs(ref[,2]-query[2])
        mtch <- which((mzdiff <= mz.tol) & (rtdiff <= rt.tol))
        # gives a unique match
        if(length(mtch)>1) mtch <- mtch[which.min(mzdiff[mtch])]
        mtch
      }
      
      mt <- apply(pqtn, 1, match.peaks, pcurve)
      
      for(i in 1:length(mt)) {
        if(length(mt[[i]])) {
          p <- qtn[,-1][,i] >= as.numeric(curve[mt[[i]], "MinInt"]) & qtn[,-1][,i] < as.numeric(curve[mt[[i]], "MaxInt"])
          if(sum(!p)) qtn[,-1][!p,i] <-  NaN
          # is that a unique match?
          #if(sum(p)) qtn[,-1][p,i] <- 10^(as.numeric(curve[mt[[i]], "Slope"])*log10(qtn[,-1][p,i]) + as.numeric(curve[mt[[i]], "Intercept"]))
          if(sum(p)) qtn[,-1][p,i] <- 10^predict.lm(modelList[[as.character(curve[mt[[i]], "Feat_pos"])]], data.frame(x=log10(qtn[,-1][p,i])))
        }
      }
      
    }
      
    if(sum(unlist(lapply(mt, length))!=0)) {
      qtn2 <- cbind(as.character(qtn[,1]), qtn[,-1][,unlist(lapply(mt, length))!=0])
      colnames(qtn2)[1] <- colnames(qtn)[1] 
      qtn2[,1] <- as.character(qtn2[,1])
      
      # adding compound names
      qtn2 <- rbind("", qtn2)
      qtn2[,-1][1, ] <- as.character(curve[unlist(mt[unlist(lapply(mt, length))!=0]),"Name"])
      
      write.csv(qtn2, file, row.names = FALSE)
    } else {
      write.csv("No feature match found.", file, row.names = FALSE)
    }
    }
  )

  EmailOutput2 <- reactive({
      input$goMail2 
    
      file = paste("quantification", '.csv', sep='') 
       name2fet <- function(fet) {
         strf <- strsplit(fet, "\\.") 
         mz <- paste(strf[[1]][1:2], collapse = ".")
         mz <- as.numeric(sub("^X", "", mz))
         rt <- as.numeric(paste(strf[[1]][3], collapse = "."))
         c(mz, rt)
      }
 
     
    if(length(input$file3)) {
      inFile <- input$file3
      #curve <- read.csv(inFile$datapath)  
      curve <- read.csv(unz(inFile$datapath, "curve_output/calibration_curve_table.csv"))
      load(unz(inFile$datapath, "curve_output/modelList.rda"))
      
      pcurve <- t(sapply(curve[,2], function(x) strsplit(as.character(x), "_")[[1]]))
      pcurve <- apply(pcurve, 2,  as.numeric)
      
      inFile <- input$file4
      qtn <- read.csv(inFile$datapath)  
      
       pqtn <- as.vector(colnames(qtn)[-1])
       pqtn <- do.call(rbind, lapply(pqtn, name2fet))
      
      mzw <- input$mz
      rtw <- input$rt
      
      match.peaks <- function(query, ref, mz.tol = mzw, rt.tol = rtw) {
        #mzdiff <- abs(ref[,1]-query[1])
        mzdiff <- (abs(ref[,1]-query[1])/query[1])*10^6
        rtdiff <- abs(ref[,2]-query[2])
        mtch <- which((mzdiff <= mz.tol) & (rtdiff <= rt.tol))
        # gives a unique match
        if(length(mtch)>1) mtch <- mtch[which.min(mzdiff[mtch])]
        mtch
      }
      
      mt <- apply(pqtn, 1, match.peaks, pcurve)
      
      for(i in 1:length(mt)) {
        if(length(mt[[i]])) {
          p <- qtn[,-1][,i] >= as.numeric(curve[mt[[i]], "MinInt"]) & qtn[,-1][,i] < as.numeric(curve[mt[[i]], "MaxInt"])
          if(sum(!p)) qtn[,-1][!p,i] <-  NaN
          # is that a unique match?
          #if(sum(p)) qtn[,-1][p,i] <- 10^(as.numeric(curve[mt[[i]], "Slope"])*log10(qtn[,-1][p,i]) + as.numeric(curve[mt[[i]], "Intercept"]))
          if(sum(p)) qtn[,-1][p,i] <- 10^predict.lm(modelList[[as.character(curve[mt[[i]], "Feat_pos"])]], data.frame(x=log10(qtn[,-1][p,i])))
        }
      }
      
    }
      
    if(sum(unlist(lapply(mt, length))!=0)) {
      qtn2 <- cbind(as.character(qtn[,1]), qtn[,-1][,unlist(lapply(mt, length))!=0])
      colnames(qtn2)[1] <- colnames(qtn)[1] 
      qtn2[,1] <- as.character(qtn2[,1])
      
      # adding compound names
      qtn2 <- rbind("", qtn2)
      qtn2[,-1][1, ] <- as.character(curve[unlist(mt[unlist(lapply(mt, length))!=0]),"Name"])
      
      write.csv(qtn2, file, row.names = FALSE)
    } else {
      write.csv("No feature match found.", file, row.names = FALSE)
    }     
          
      from = "<mysterious-cloud@dal.ca>"
      to = paste0("<", input$mail2, ">")
      
      control=list(smtpServer="ASPMX.L.GOOGLE.COM")
      subject <- "Quantification table"
      body <- "Here are your calibration curves, please let us know if the formatting is correct."
      #sendmail(from=from,to=to,subject=subject,msg=body,control=control)
      attachmentPath <- file 
      attachmentName <- file 
      attachmentObject <- mime_part(x=attachmentPath,name=attachmentName)
      bodyWithAttachment <- list(body,attachmentObject)
      
      sendmail(from=from,to=to,subject=subject,msg=bodyWithAttachment,control=control)
  })
  
  observeEvent(input$goMail2, {
    EmailOutput2()
  })
 
   
  
})
