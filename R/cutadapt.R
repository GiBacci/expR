.Cutadapt <- setClass("Cutadapt",
                     slots = c(
                       adapter = "character",
                       output = "character",
                       input = "character"
                     ))

validCutadaptObject <- function(object){
  if(object@output == object@input){
    paste("Input and output are the same!")
  }else{
    TRUE
  }
}

setValidity("Cutadapt", validCutadaptObject)

cutadapt <- function(adapter, output = "output.txt", input = "input.txt"){
  .Cutadapt(adapter = adapter, output = output, input = input)
}

c <- cutadapt("ACTGC")
