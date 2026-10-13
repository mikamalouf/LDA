library(KernSmooth)
library(cgam)


get_fit<-function(std,conc,type){
    max <- max(std[which(is.finite(std))])

    # if (is.na(std[6])==TRUE || is.na(std[5])==TRUE) {

    # } else if (std[6] < std[5]) {
        std <- std[-6]
        conc <- conc[-6]
    # }
        # std <- std[-1]
        # conc <- conc[-1]
    print(std)
    print(conc)

    len_std <- length(std)
    fit.data <- data.frame("MFI"=std,"conc"=log(conc))
    
    if (type=="cgam"){
        fit<-cgam(MFI ~ 1+s.incr(conc,numknots=4), data=fit.data, family=gaussian)
        y1 <-predict.cgam(fit,newData = data.frame(conc=x1))$fit
    } else if (type=="kernsmooth"){
        fit <- locpoly(
                y = fit.data$MFI,
                x = fit.data$conc,
                degree = 1,
                bandwidth = 1.2,
                kernel = "Epanechnikov",
                gridsize = 100
            )
        x1 <- fit$x
        y1 <- fit$y
    } else if (type=="poly"){
        fit<-lm(MFI ~ poly(conc, 3, raw = TRUE), data = fit.data)
        params <- as.list(coef(fit))
        scale <- params$Scale
        mid <- params$Xmid
        max <- params$maxMFI
        x1 <- seq(min(fit.data$conc),max(fit.data$conc),length=100)
        y1 <- predict(fit, newdata=list(conc=x1))
    } else if (type=="nls"){
        cc1 <- try( starting.values <- getInitial(MFI ~ SSlogis(conc,maxMFI,Xmid,Scale),data=fit.data), silent=T)
        cc2 <- try( nls(MFI ~ (maxMFI / (1 + exp((Xmid - conc)/Scale))),data=fit.data, start=starting.values), silent=T)
        if (is(cc1,"try-error")|is(cc2,"try-error")) {
        if (is.finite(std[len_std])==T) std_last <- std[len_std] else std_last <- std[len_std-1]

        starting.values <- expand.grid(maxMFI=seq(0, std_last*2.5, len = 20),
                                        Xmid=seq(-50, 0, len = 20),
                                        Scale=seq(0, 50, len = 40))

        fit <- nls2(MFI ~ (maxMFI / (1 + exp((Xmid - conc)/Scale))), algorithm="brute-force",data=fit.data, start=starting.values)

        } else {

        starting.values <- getInitial(MFI ~ SSlogis(conc,maxMFI,Xmid,Scale),data=fit.data)
        fit <- nls(MFI ~ (maxMFI / (1 + exp((Xmid - conc)/Scale))),data=fit.data, start=starting.values)

        }

        params <- as.list(coef(fit))

        scale <- params$Scale
        mid <- params$Xmid
        max <- params$maxMFI

        x1 <- seq(min(fit.data$conc),max(fit.data$conc),length=100)
        y1 <- predict(fit, newdata=list(conc=x1))
    }
    list(x=x1,y=y1)
}

get_conc<-function(dilutions){
    sapply(strsplit(dilutions, split = "/"),
                        function(x) as.numeric(x[1]) / as.numeric(x[2]))
}

normalise_plates <- function(all_plates,dilutions,ag_list,fit_type,progress){
    fitted.data<-list()
    plots<-list()
    for (k in 1:length(all_plates)) {

        plate_lab <- names(all_plates)[k]
        
        standards <- get.standard(all_plates[[k]], std_label = "(?i)CP3|Std Curve", dilutions, n_points=length(dilutions))
        
        fit.mat <- as.data.frame(matrix(ncol=length(ag_list)+1,nrow=100))
        fit.mat[,1] <- seq(1,100,1)
        colnames(fit.mat) <- c("conc_no",ag_list)
        
        for (i in 1:length(ag_list) ) {
            ag<-ag_list[i]
            print(ag)
            std <- standards[,ag_list[i]]
            conc <- get_conc(dilutions)

            len_std <- length(std)
            
            fit<-get_fit(std,conc,type=fit_type)
            x1<-fit$x
            y1<-fit$y
            fit.mat[,i+1] <- y1


            plot(log(conc),std, type="p", pch=19, col="dodgerblue4",ylab="MFI", xlab="Dilution factor",axes=F,ylim=c(0,max(y1)),
                    xlim=c(min(log(conc)),max(log(conc))+1),main=paste0(ag_list[i],"Standard curve\n",plate_lab), cex.main=0.8, cex.lab=0.8 )
                axis(1, at=c(log(conc),max(log(conc))+1) , labels=c(dilutions[1:len_std],""), cex.axis=0.8)
                axis(2, cex.axis=0.8)
                lines(x1,y1,type="l",lwd=2, col="dodgerblue4")

            plots[[plate_lab]][[ag]] <- recordPlot()
        }
        fitted.data[[k]]<-fit.mat
        if (! is.null(progress)){
            progress$set(value = k)
        }
    }
    names(fitted.data)<-names(all_plates)
    list(data=fitted.data,plots=plots)
}

loess_adjustment<-function(all_plates,fitted_data,ref_plate,ag_list,progress){
    fit.mat_clean<-fitted_data
    ref_fit <- list(fit.mat_clean[[ref_plate]])
    loess_plate <- list()
    plots <-  list()

    for (i in 1:length(fit.mat_clean)){
        plate_lab <- names(fit.mat_clean)[[i]]
        plots[[plate_lab]]<-list()
        print(plate_lab)
        plate_fit <- fit.mat_clean[[plate_lab]]

        h1_raw <- all_plates[[plate_lab]] #[-grep("CP3",all_plates[[plate_lab]][,"Sample"]),]
        h1_raw <-  h1_raw[-grep("Background0", h1_raw$Sample),]
        idx<- grep("NEG", h1_raw$Sample)
        if (length(idx)>0){
        h1_raw <-  h1_raw[-grep("NEG", h1_raw$Sample),]
        }
        idx<- grep("WHO", h1_raw$Sample)
        if (length(idx)>0){
        h1_raw<-   h1_raw[-grep("WHO", h1_raw$Sample),]
        }
        
        
        #If plate has this control
        #h1_raw <-  h1_raw[-grep("N.am", h1_raw$Sample),]
        loess_plate[[i]] <- norm.loess(
            data = h1_raw,
            ref_fit = ref_fit,
            plate_fit = plate_fit,
            ag_list = ag_list,
            plot_fit = FALSE
        )

        for (j in ag_list) {

            h1_data <- h1_raw[,j]
            h2_data <- as.numeric(loess_plate[[i]][,j])

            max <- max(c(h1_data,h2_data),na.rm=T)
            min <- min(c(h1_data,h2_data),na.rm=T)
            equal <- seq(0,max,by=1)

            plot(
                h1_data,
                h2_data,
                cex = 0.8,
                main = paste0(plate_lab,"\n",j),
                ylim = c(min,max*1.1),
                xlim = c(min,max*1.1), 
                xlab = "Raw data", 
                ylab = "Norm data"
            )
            lines(equal,equal,lty=3,col="red",lwd=0.5)
            plots[[plate_lab]][[j]] <- recordPlot()

        }
        if (! is.null(progress)){
            progress$set(value = i)
        }
    }
    names(loess_plate)<-names(fit.mat_clean)

    list(
        data = loess_plate,
        plots = plots
    )
}


plot.ag.plates.raw <- function(data, ag, dilutions) {
    conc<-get_conc(dilutions)
    tmp <- as.data.frame(sapply(names(data),function(p){
        as.numeric(get.standard(data[[p]], std_label = "(?i)CP3|Std Curve", dilutions, n_points=length(dilutions))[[ag]][1:5])
    }))
    conc<-conc[1:nrow(tmp)]
    print(tmp)
    plot(1:nrow(tmp),tmp$MP01,ylim=c(0,max(tmp)),type="l",main=ag,xaxt="n",ylab="MFI",xlab="Dilution")
    axis(side=1,at=1:nrow(tmp),labels=dilutions[1:nrow(tmp)])
    apply(tmp,2,function(x){
        # print(x)
        lines(1:nrow(tmp),x)
    })
}

plot.ag.plates.norm <- function(data, ag, dilutions) {
    conc<-get_conc(dilutions)
    tmp <- as.data.frame(sapply(names(data),function(p){
        data[[p]][ag]
    }))
    plot(1:100,tmp$MP01,ylim=c(0,max(tmp)),type="l")
    apply(tmp,2,function(x){
        lines(1:100,x)
    })
}