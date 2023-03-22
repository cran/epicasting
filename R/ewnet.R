#' Ewnet: An Ensemble Wavelet Neural Network for Forecasting and Epicasting
#'
#' @param ts A numeric vector or time series
#' @param Waveletlevels An integer specifying the levels of decomposition. The default
#' is set to floor(log(length(ts))).
#' @param MaxARParam An integer indicating the maximum lagged observations to be included
#' in the neural network. The default is selected based on AIC using linear AR process.
#' @param boundary  A character string indicating which boundary method to use.
#' boundary = "periodic" (default) and boundary = "reflection".
#' @param FastFlag A logical flag which, if true (default), indicates that the pyramid
#' algorithm is computed with an internal C function. Otherwise, only R code is used
#' in all computations.
#' @param NForecast An integer specifying the forecast horizon.
#' @param NVal An integer indicating the size of validation set. Default is set to 0.
#' @param measure The performance metric used for selecting the best value of \code{MaxARParam}
#' based on validation set. Defaults to Metrics::mase.
#' @param PI A logical flag which, if true generates the confidence interval for the
#' forecast horizon. Default is set to false.
#' @param xreg_train Optionally, a vector or matrix of external regressors, which
#' must have the same number of rows as \code{ts}. Must be numeric.
#' @param xreg_test Optionally, a vector or matrix of external regressors, which
#' must have the same number of rows as \code{NForecast} to be used for the forecast.
#' Must be numeric.
#' @param ret_fit A logical flag specifying that the fitted values of the model on the
#' training set should be returned if true, otherwise, false (default).
#' @import wavelets stats
#' @importFrom Metrics mase
#' @importFrom forecast nnetar
#' @importFrom forecast forecast
#'
#' @return The parameters of the fitted model indicating the number of lagged observations
#' included in the model and the number of nodes in the hidden layer. The forecast of the
#' time series of size \code{NForecast} is generated along with the optional output of
#' fitted values (\code{ret_fit} = TRUE) and confidence interval (\code{PI} = TRUE) for the forecast.
#'
#' @export
#'
#' @author Madhurima Panja and Tanujit Chakraborty
#'
#' @references \itemize{
#' \item Panja, M., Chakraborty, T., Kumar, U., & Liu, N. (2022).
#'       Epicasting: An ensemble wavelet neural network (ewnet) for forecasting epidemics.
#'       arXiv preprint arXiv:2206.10696. \url{https://arxiv.org/abs/2206.10696}
#'
#' \item Panja, M., Chakraborty, T., Nadim, S. S., Ghosh, I., Kumar, U., & Liu, N. (2023).
#'       An ensemble neural network approach to forecast Dengue outbreak based on climatic condition.
#'       Chaos, Solitons & Fractals, 167, 113124.}
#'
#' @examples
#' ewnet(ts = datasets::lynx, MaxARParam = 1, NForecast = 3)

ewnet <- function(ts,Waveletlevels = floor(log(length(ts))),MaxARParam,boundary = "periodic",FastFlag = TRUE,NForecast,NVal = 0, measure = Metrics::mase, PI = FALSE, xreg_train = NULL, xreg_test = NULL, ret_fit = FALSE){
  WaveletFittingnar<- function(ts,Waveletlevels,MaxARParam,boundary,FastFlag,NForecast, xreg_train = NULL, xreg_test = NULL)
  {
    WaveletFitting <- function(ts,Wvlevels,bndry,FFlag)
    {
      mraout <- wavelets::modwt(ts, filter='haar', n.levels=Wvlevels,boundary=bndry, fast=FFlag)
      WaveletSeries <- cbind(do.call(cbind,mraout@W),mraout@V[[Wvlevels]])
      return(list(WaveletSeries=WaveletSeries,WVSeries=mraout))
    }
    WS <- WaveletFitting(ts=ts,Wvlevels=Waveletlevels,bndry=boundary,FFlag=FastFlag)$WaveletSeries
    AllWaveletForecast <- NULL;AllWaveletPrediction <- NULL

    for(WVLevel in 1:ncol(WS))
    {
      ts <- NULL
      ts <- WS[,WVLevel]
      WaveletNARFit <- forecast::nnetar(y=stats::as.ts(ts), p = MaxARParam, repeats = 500, xreg = xreg_train)
      lags = WaveletNARFit$p ## New
      WaveletNARPredict <- WaveletNARFit$fitted
      WaveletNARForecast <- forecast::forecast(WaveletNARFit, h=NForecast, xreg = xreg_test)
      AllWaveletPrediction <- cbind(AllWaveletPrediction,WaveletNARPredict)
      AllWaveletForecast <- cbind(AllWaveletForecast,as.matrix(WaveletNARForecast$mean))
    }
    Finalforecast <- rowSums(AllWaveletForecast,na.rm = T)
    FinalPrediction <- rowSums(AllWaveletPrediction,na.rm = T)
    return(list(Finalforecast=Finalforecast,FinalPrediction=FinalPrediction, lags = lags))
  }
  n_test = NForecast
  if(NVal == 0){
    fit_ewnet = WaveletFittingnar(ts(ts), Waveletlevels = floor(log(length(ts))), boundary = "periodic",
                                  FastFlag = TRUE, MaxARParam, NForecast, xreg_train = NULL, xreg_test = NULL)
    fit_ewnet$Finalforecast[fit_ewnet$Finalforecast<0] = 0
  }else{
    train_val = subset(ts(ts),  end= length(ts(ts))-NVal)
    val = subset(ts(ts),  start= length(ts)-NVal+1)
    n_val = length(val)
    model_smry <- data.frame()
    if(is.null(xreg_train)){
      for(p in 1:MaxARParam){
        fit_ewnet_val = WaveletFittingnar(ts(train_val), Waveletlevels = floor(log(length(train_val))), boundary = "periodic",
                                          FastFlag = TRUE, MaxARParam = p, NForecast = n_val, xreg_train = NULL, xreg_test = NULL)
        fore_ewnet_val = as.data.frame(fit_ewnet_val$Finalforecast, h = n_val)
        ewnet_val_score = measure(val, fore_ewnet_val$`fit_ewnet_val$Finalforecast`)
        ewnet_val_evaluation = data.frame(score = ewnet_val_score, p)
        model_smry <- rbind(model_smry, ewnet_val_evaluation)
      }
    }else{
      xreg_train_val = subset(ts(xreg_train),  end= length(ts(xreg_train))-NVal)
      xreg_val = subset(ts(xreg_train),  start= length(xreg_train)-NVal+1)
      WaveletFittingval<- function(ts,Waveletlevels,MaxARParam,boundary,FastFlag,NForecast, xreg_train_val = NULL, xreg_val = NULL)
      {
        WaveletFitting <- function(ts,Wvlevels,bndry,FFlag)
        {
          mraout <- wavelets::modwt(ts, filter='haar', n.levels=Wvlevels,boundary=bndry, fast=FFlag)
          WaveletSeries <- cbind(do.call(cbind,mraout@W),mraout@V[[Wvlevels]])
          return(list(WaveletSeries=WaveletSeries,WVSeries=mraout))
        }
        WS <- WaveletFitting(ts=ts,Wvlevels=Waveletlevels,bndry=boundary,FFlag=FastFlag)$WaveletSeries
        AllWaveletForecast <- NULL;AllWaveletPrediction <- NULL

        for(WVLevel in 1:ncol(WS))
        {
          ts <- NULL
          ts <- WS[,WVLevel]
          WaveletNARFit <- forecast::nnetar(y=stats::as.ts(ts), p = MaxARParam, repeats = 500, xreg = xreg_train_val)
          lags = WaveletNARFit$p
          WaveletNARPredict <- WaveletNARFit$fitted
          WaveletNARForecast <- forecast::forecast(WaveletNARFit, h=NForecast, xreg = xreg_val)
          AllWaveletPrediction <- cbind(AllWaveletPrediction,WaveletNARPredict)
          AllWaveletForecast <- cbind(AllWaveletForecast,as.matrix(WaveletNARForecast$mean))
        }
        Finalforecast <- rowSums(AllWaveletForecast,na.rm = T)
        FinalPrediction <- rowSums(AllWaveletPrediction,na.rm = T)
        return(list(Finalforecast=Finalforecast,FinalPrediction=FinalPrediction, lags = lags))
      }
      for(p in 1:MaxARParam){
        fit_ewnet_val = WaveletFittingval(ts(train_val), Waveletlevels = floor(log(length(train_val))), boundary = "periodic",
                                          FastFlag = TRUE, MaxARParam = p, NForecast = n_val, xreg_train_val = NULL, xreg_val = NULL)
        fore_ewnet_val = as.data.frame(fit_ewnet_val$Finalforecast, h = n_val)
        ewnet_val_score = measure(val, fore_ewnet_val$`fit_ewnet_val$Finalforecast`)
        ewnet_val_evaluation = data.frame(score = ewnet_val_score, p)
        model_smry <- rbind(model_smry, ewnet_val_evaluation)
      }
    }
    final = model_smry[which.min(model_smry$score),]
    fit_ewnet = WaveletFittingnar(ts(ts), Waveletlevels = floor(log(length(ts))), boundary = "periodic",
                                  FastFlag = TRUE, MaxARParam = final$p, NForecast = n_test, xreg_train = NULL, xreg_test = NULL)
    fit_ewnet$Finalforecast[fit_ewnet$Finalforecast<0] = 0
  }
  if (isTRUE(PI)){
    upper = fit_ewnet$Finalforecast + 1.5*stats::sd(fit_ewnet$Finalforecast)
    lower = fit_ewnet$Finalforecast - 1.5*stats::sd(fit_ewnet$Finalforecast)
    forecast = list("Parameters" = c(fit_ewnet$lags, round((fit_ewnet$lags+1)/2)),
                    "Fitted" = fit_ewnet$FinalPrediction,
                    "Forecast" = fit_ewnet$Finalforecast,
                    "Lower Interval" = lower,
                    "Upper Interval" = upper)
  }else{
   forecast = list("Parameters" = c(fit_ewnet$lags, round((fit_ewnet$lags+1)/2)),
                    "Fitted" = fit_ewnet$FinalPrediction,
                    "Forecast" = fit_ewnet$Finalforecast)
  }
  if(isTRUE(ret_fit)){
    forecast = forecast
  }else{
    forecast = forecast[-2]
  }
  return(forecast)
}
