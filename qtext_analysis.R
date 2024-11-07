############## Code for conducting analysis ############################ ############# ############
#### "WHAT PATIENTS CALL GYNECOLOGICAL CONDITIONS: A QUALITATIVE STUDY". November, 2024 ##########
################################################################################################## 

library(udpipe)
library(tidyverse)
library(tidytext) 
library(ngram)
library(spacyr)
library(word2vec)
library(stringdist)
library(igraph)
library(ggraph)

set.seed(85566)
stopwords <- tidytext::stop_words$word


theme_tuf2 <- ggthemes::theme_tufte() +
  theme(
    axis.line = element_line(color = 'black'),
    axis.title.x = element_text(vjust = -0.3), 
    axis.title.y = element_text(vjust = 0.8),
    legend.background = element_blank(), 
    legend.key = element_blank(), 
    legend.title = element_blank(),
    panel.background = element_blank(), 
    panel.border = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_blank(),
    strip.background = element_blank()
  )

# Aux functions:

get_vocs <- function(x){
  
  r1 <- spacy_parse(x, tag = TRUE, pos = TRUE)
  nouns <- subset(r1, pos == 'NOUN')
  nouns <- split(nouns$token, nouns$doc_id)$text1
  r2 <- spacy_parse(x, tag = TRUE, pos = TRUE)
  adjs <- subset(r2, pos == 'ADJ')
  adjs <- split(adjs$token, adjs$doc_id)$text1
  r3 <- spacy_parse(x, tag = TRUE, pos = TRUE)
  advs <- subset(r2, pos == 'ADV')
  advs <- split(advs$token, advs$doc_id)$text1
  
  return(c(nouns, adjs, advs))
}

compare_grpterm <- function(df, term, grp){
  
  dfg <- df %>%
    select(c(grp, "Q1_voclength", "Q1_vocs", "Q1_clean")) %>%
    group_by_at(c(grp)) %>% 
    mutate(grptxt = paste0(Q1_vocs, collapse = " ")) %>% 
    mutate(termfreq = str_count(term, grptxt))
  
  dfg2 <- df %>%
    group_by_at(c(grp)) %>% 
    summarise(tot_textlen = sum(unlist(Q1_voclength)))
  
  if(grp=="Race"){dfg$grp <- as.factor(dfg$Race)}
  if(grp=="Age"){dfg$grp <- as.factor(dfg$Age)}
  if(grp=="Num_childbirths"){dfg$grp <- as.factor(dfg$Num_childbirths)}
  if(grp=="Job_field"){dfg$grp <- as.factor(dfg$Job_field)}
  
  dfg2$grptxt[1] <- filter(dfg, grp==levels(dfg$grp)[1])$grptxt[1]
  dfg2$grptxt[2] <- filter(dfg, grp==levels(dfg$grp)[2])$grptxt[1]
  
  dfg2$termfreq[1] = str_count(dfg2$grptxt[1], term)/dfg2$tot_textlen[1]
  dfg2$termfreq[2] = str_count(dfg2$grptxt[2], term)/dfg2$tot_textlen[2]
  
  dfg2$term <- term
  
  return(dfg2)
  
}

check_term <- function(df, w, grp, model, simthresh=0.85){
  
  dft <- df %>%
    
    select(c(grp, "Q1_clean")) 
  
  dft$Term <- w
  dft$Mentioned <- ifelse(str_detect(dft$Q1_clean, w), "Yes","No")
  
  for (i in nrow(dft)){
    if(dft$Mentioned[i] == "No"){
      swd=unique(filter(dfsynons, term==w)$synons)
      if(length(swd) > 0) {
        for(sw in swd){
          if(str_detect(dft$Q1_clean[i], sw))
          {dft$Mentioned[i] <- "Yes"
          print("was there alt1")
          }
        }
      }
    }
  }
  
  for (i in nrow(dft)){
    if(dft$Mentioned[i] == "No"){
      txtwords <- str_split(dft$Q1_clean[i], " ")[[1]]
      txtwords <-  txtwords[nchar(txtwords) > 1]
      for(tw in txtwords){
        if(stringdist::stringsim(w, tw,
                                 method="jw") >= simthresh){
          {dft$Mentioned[i] <- "Yes"
          print("was there alt2")
          }
        }
      }
    }
  }
  
  
  
  for (i in nrow(dft)){
    if(dft$Mentioned[i] == "No"){
      if(w %in% unique(dfsynons$synons)){
        swf <- filter(dfsynons, synons==w)$term[1]
        if(str_detect(dft$Q1_clean[i], swf))
        {dft$Mentioned[i] <- "Yes"
        print("was there alt3")
        }
      }
    }
  }
  
  for (i in nrow(dft)){
    if(dft$Mentioned[i] == "No"){
      
      embedding <- predict(model, c(w), type = "embedding")
      if(!is.na(embedding)[,1]){
        lookslike <- predict(model, c(w), type = "nearest", top_n=1000)
        lookslike <- as.data.frame(lookslike[w])
        colnames(lookslike) <- c("term1", "term2", "similarity", "rank")
        for(sw in lookslike$term2){
          if(str_detect(dft$Q1_clean[i], sw))
          {dft$Mentioned[i] <- "Yes"
          print("OH 4")
          print(w)
          print(sw)
          print("was there alt4")
          }
        }
      }
    }
  }
  
  
  return(select(dft, c(grp, Term, Mentioned)))
  
}

test_term <- function(dft, clinterms, grp, model, Qcol){
  
  gtlist = list()
  
  for(i in 1:length(clinterms)){
    term <- clinterms[[i]]
    dfterm <- filter(dft, Term==term)
    dfterm <- select(dfterm, c(grp, "Mentioned"))
    gtlist[[i]] <- gtsummary::tbl_summary(dfterm,
                                          by = as.name(grp),
                                          percent = "column", 
                                          statistic = list(all_continuous() ~ "{median} ({p25}, {p75})",
                                                           all_categorical() ~ "{n} ({p}%)"),
                                          label = list(Mentioned ~ " ")
                                          
    ) %>%
      gtsummary::modify_header(label = " ") %>%
      gtsummary::bold_labels() %>%
      gtsummary::modify_header(all_stat_cols() ~ "**{level}**, N = {n}") %>%
      gtsummary::add_p()
  }
  
  
  gt <- gtsummary::tbl_stack(tbls = gtlist, group_header = clinterms) 
  gt <- as_gt(gt) %>%
    gt::tab_options(table.font.names = "Times New Roman",
                    table.font.size = 14) %>%
    gt::tab_header(title = paste0(str_to_title(str_replace(Qcol, "Qs", "")),
                                  ": Mentioned Terms by: ",
                                  str_replace(str_to_title(grp), "_", " ")))
  # }
  return(gt)
}

prep_Qcol <- function(df, Qcol){
  
  df$Q1 <- df[,Qcol]
  
  puncts <- tm::removePunctuation(df$Q1)
  
  
  for (sw in stopwords){
    sw <- paste0(" ", sw, " ")
    df$Q1_clean <- str_replace_all(df$Q1, sw, " ")
  }
  
  
  df$Q1_clean <- tolower(df$Q1_clean)
  df$Q1_clean <- str_replace_all(df$Q1_clean, " i ", " ")
  df$Q1_clean <- str_replace_all(df$Q1_clean, " a ", " ")
  df$Q1_clean <- str_replace_all(df$Q1_clean, "the ", " ")
  df$Q1_clean <- str_replace_all(df$Q1_clean, " of", " ")
  df$Q1_clean <- str_replace_all(df$Q1_clean, " and ", " ")
  df$Q1_clean <- str_replace_all(df$Q1_clean, " they ", " ")
  df$Q1_clean <- str_replace_all(df$Q1_clean, " to ", " ")
  df$Q1_clean <- str_replace_all(df$Q1_clean, " you ", " ")
  df$Q1_clean <- str_replace_all(df$Q1_clean, " your ", " ")
  df$Q1_clean <- str_replace_all(df$Q1_clean, " but ", " ")
  df$Q1_clean <- str_replace_all(df$Q1_clean, " as ", " ")
  df$Q1_clean <- str_replace_all(df$Q1_clean, " like ", " ")
  df$Q1_clean <- str_replace_all(df$Q1_clean, " at ", " ")
  df$Q1_clean <- str_replace_all(df$Q1_clean, " or ", " ")
  df$Q1_clean <- str_replace_all(df$Q1_clean, "so ", " ")
  df$Q1_clean <- str_replace_all(df$Q1_clean, " for ", " ")
  df$Q1_clean <- str_replace_all(df$Q1_clean, " with ", " ")
  df$Q1_clean <- str_replace_all(df$Q1_clean, " or ", " ")
  df$Q1_clean <- str_replace_all(df$Q1_clean, "yeah ", " ")
  df$Q1_clean <- str_replace_all(df$Q1_clean, "um ", " ")
  df$Q1_clean <- str_replace_all(df$Q1_clean, " in ", " ")
  df$Q1_clean <- str_replace_all(df$Q1_clean, " my ", " ")
  df$Q1_clean <- str_replace_all(df$Q1_clean, " that ", " ")
  df$Q1_clean <- str_replace_all(df$Q1_clean, " have ", " ")
  df$Q1_clean <- str_replace_all(df$Q1_clean, " like ", " ")
  df$Q1_clean <- str_replace_all(df$Q1_clean, " know ", " ")
  df$Q1_clean <- str_replace_all(df$Q1_clean, " think ", " ")
  df$Q1_clean <- str_replace_all(df$Q1_clean, " it ", " ")
  df$Q1_clean <- str_replace_all(df$Q1_clean, " about ", " ")
  df$Q1_clean <- str_replace_all(df$Q1_clean, " just ", " ")
  df$Q1_clean <- str_replace_all(df$Q1_clean, " then ", " ")
  df$Q1_clean <- str_replace_all(df$Q1_clean, " more ", " ")
  
  
  
  
  
  
  df$Q1_clean <- tm::removePunctuation(df$Q1_clean)
  df$Q1_clean <- str_replace_all(df$Q1_clean, "  ", " ")
  
  df$Q1_vocs <- lapply(df$Q1_clean, function(x) {unlist(get_vocs(x))})
  df$Q1_vocs <- lapply(df$Q1_vocs, function(x) glue::glue_collapse(x, sep= " "))
  
  df$Q1_voclength <- lapply(df$Q1_vocs, function(x) {length(str_split(x, " ")[[1]])})
  
  return(df)
}

get_word_freq <- function(dfW, txtfilename="", grp="all"){
  
  dfW2 <- apply(select(dfW, "Q1_vocs"), 1, paste, collapse="")
  
  if(grp != "all"){
    
    dfW2 <- as.data.frame(dfW2)
    
    dfW2$grp <- dfW[grp]
    
    colnames(dfW2) <- c("txt", "grp")
    
    dfWG <- dfW2 %>% 
      group_by(grp) %>% 
      mutate(grptxt = paste(txt, collapse = " ")
      ) 
    
    dfWG  <- dfWG [!duplicated(dfWG [,c('grptxt')]),][,-1]
    colnames(dfWG) <- c("grp", "txt")
    freq1 <- table(str_split(dfWG$txt[1], " "))
    freq1 <- sort(freq1, decreasing=TRUE)
    freq1 <- as.data.frame(freq1)
    write_csv(freq1,
              paste("results/", txtfilename, "_", unique(dfWG$grp)[[1]][1], ".csv" ))
    
    freq2 <- table(str_split(dfWG$txt[2], " "))
    freq2 <- sort(freq2, decreasing=TRUE)
    freq2 <- as.data.frame(freq2)
    write_csv(freq2,
              paste("results/", txtfilename, "_", unique(dfWG$grp)[[1]][2], ".csv" ))
    
    freq <- list(freq1, freq2)

  }
  
  if(grp=="all"){
    
    dfW2 <- paste(dfW2 , collapse = " ")
    freq <- table(str_split(dfW2, " "))
    freq <- sort(freq, decreasing=TRUE)
    freq <- as.data.frame(freq)
    write_csv(freq, 
              paste("results/", txtfilename, "_", "all", ".csv" ))
  }
  
  return(freq)
}

get_sents <- function(){
  
  s1<- syuzhet::get_nrc_sentiment(as.vector(tAll$Var1))
  s1 <- as.data.frame(colMeans(s1))
  s1$sentiment <- rownames(s1)
  s1$var <- "All"
  colnames(s1) <- c("score", "sentiment", "var")
  
  s2 <- syuzhet::get_nrc_sentiment(as.vector(tPer$Var1))
  s2 <- as.data.frame(colMeans(s2))
  s2$sentiment <- rownames(s2)
  s2$var <- "Period"
  colnames(s2) <- c("score", "sentiment", "var")
  
  s3 <- syuzhet::get_nrc_sentiment(as.vector(tVag$Var1))
  s3<- as.data.frame(colMeans(s3))
  s3$sentiment <- rownames(s3)
  s3$var <- "Vaginal"
  colnames(s3) <- c("score", "sentiment", "var")
  
  s4 <- syuzhet::get_nrc_sentiment(as.vector(tUri$Var1))
  s4 <- as.data.frame(colMeans(s4))
  s4$sentiment <- rownames(s4)
  s4$var <- "Urine"
  colnames(s4) <- c("score", "sentiment", "var")
  
  s5 <- syuzhet::get_nrc_sentiment(as.vector(tBow$Var1))
  s5 <- as.data.frame(colMeans(s5))
  s5$sentiment <- rownames(s5)
  s5 $var <- "Bowel"
  colnames(s5) <- c("score", "sentiment", "var")
  
  s6 <- syuzhet::get_nrc_sentiment(as.vector(tProl$Var1))
  s6 <- as.data.frame(colMeans(s6))
  s6$sentiment <- rownames(s6)
  s6$var <- "Prolapse"
  colnames(s6) <- c("score", "sentiment", "var")
  
  s7 <- syuzhet::get_nrc_sentiment(as.vector(tFibr$Var1))
  s7 <- as.data.frame(colMeans(s7))
  s7$sentiment <- rownames(s7)
  s7$var <- "Fibrosis"
  colnames(s7) <- c("score", "sentiment", "var")
  
  s <- do.call(rbind, list(s1, s2, s3, s4, s5, s6, s7))
  s$sentiment <- str_to_title(s$sentiment)
  
  return(s)
}


######### APPLYING FUNCTIONS

## Read data
df <- read_csv(file_name)


# Define terms of clinical interest
clintermsPER <- c("pain", "cramp" ,"heat", "cancer", "fatigue", "bleeding", 
               "mood", "headache", "hair", "incontinence",
               "back", "work", "menopause",
               "bloating", "weight", "hormone",
               "endometriosis", "imbalance")

clintermsURIN <- c("pain", "bleeding", "incontinence", "urine", "leakage", 
                   "leak", "cough", "sneeze",
                   "incontient", "urinary", "pee", "bladder")

clintermsBOW <- c("pain", "bleeding", "stool", "poop", "ibs", "constipation", "hemorrhoid",
                  "diarrhea","constipated", "leakage", "incontinence", "incontinent")

clintermsVAG <- c("pain", "bleeding", "vagina", "yeast", "dry", "dryness", "discharge", "smell",
                  "infection", "itch", "cervix")

clintermsPROL <- c("pain", "bleeding", "uterus", "prolapse",
                   "birth", "protruding", "fall", "bulge")

clintermsFIBR <- c("pain", "bleeding", "fibroid", "cysts",
                   "cancer", "tumor", "polyp")

clintermsMENOP <- c("menopause", clintermsPER)

###### Medical Word2Vec to find semantically similar words:

loadWord2Vec <- function(build=FALSE,
                         filepath="wordict/word2vecMed.RDS"){
  
  if(build==TRUE){
  xstops <- c("in", "on", "over", "under",
              "is", "are", "were", "was",
              "a",  "had",
              "the",  "about",
              "of",  
              "and",  
              "they",  
              "to",  
              "you",  
              "your",  
              "but",  
              "as",  
              "like",  
              "at",  
              "or",  
              "so",  
              "for",  
              "with",  
              "or", 
              "yeah",  
              "um",  
              "in",  
              "my",  
              "that",  
              "have",  
              "like",  
              "know",  
              "think",  
              "it",  
              "about",  
              "just",  
              "them",
              "me")
  # https://www.kaggle.com/datasets/chaitanyakck/medical-text/data
  
  x <- readLines("wordict/train.dat")
  
  x <- tolower(x)
  x <- tm::removePunctuation(x)
  x <- tm::removeNumbers(x)
  x <- tm::removeWords(x, stopwords)
  
  model <- word2vec(x = x, type = "skip-gram", iter = 20, stopwords = xstops, min_count=3)
  
  #embedding <- as.matrix(model)
  saveRDS(model, "wordict/word2vecMed.RDS")
  } else {model <- readRDS("wordict/word2vecMed.RDS")}
  
  return(model)
}

model <- loadWord2Vec(build=TRUE)


# Topics: periods, fibroids, menopause, urinary incontinence, bowel, and vaginal

# group questions:


dfQs <- df

dfQs$periodQs <- paste(dfQs$period1, dfQs$period2, dfQs$period3,
                       dfQs$period24, dfQs$period25, dfQs$period26)

dfQs$urinQs <- paste(dfQs$urin4, dfQs$urin5, dfQs$urin13, dfQs$urin14, 
                     dfQs$urin15, dfQs$urin16
                     )

dfQs$fibrQs <- paste(dfQs$fibr22, dfQs$fibr23)
dfQs$menopQs <- paste(dfQs$menop27, dfQs$menop28)
dfQs$bowQs <- paste(dfQs$bow7, dfQs$bow8, dfQs$bow18, dfQs$bow19, dfQs$bow20,
                    dfQs$bow21)
dfQs$vagQs <- dfQs$vag6
dfQs$prolapseQs <- paste(dfQs$prol9, dfQs$prol10, dfQs$prol11, dfQs$prol12)

dfQs <- select(dfQs, c("Race", "Age", "Num_childbirths", "Job_field", 
                       colnames(dfQs)[33:39]  ))

### Synonyms dict (manual):
term <- "pain"
synons <- c("ache", "painful", "discomfort", "sore")
dfsynons <- data.frame(term, synons)

term <- "cancer"
synons <- c("tumor", "tumorous")
dfsynons <- rbind(dfsynons, data.frame(term, synons))

term <- "mood"
synons <- c("irritable", "irritability", "swings", "emotional")
dfsynons <- rbind(dfsynons, data.frame(term, synons))

term <- "tired"
synons <- c("fatigue", "exhasut", "energy", "exhaust", "sleep")
dfsynons <- rbind(dfsynons, data.frame(term, synons))

term <- "cravings"
synons <- c("food", "crave", "appetite", "craving")
dfsynons <- rbind(dfsynons, data.frame(term, synons))

term <- "weight"
synons <- c("bodyweight", "metabolism")
dfsynons <- rbind(dfsynons, data.frame(term, synons))

term <- "hormonal"
synons <- c("hormone", "imbalance", "estrogen")
dfsynons <- rbind(dfsynons, data.frame(term, synons))

term <- "work"
synons <- c("performance", "job", "productivity", "career")
dfsynons <- rbind(dfsynons, data.frame(term, synons))

term <- "nausea"
synons <- c("vomit", "dizzy", "stomach", "stomachache")
dfsynons <- rbind(dfsynons, data.frame(term, synons))

term <- "bleeding"
synons <- c("blood", "bleed", "period", "spotting")
dfsynons <- rbind(dfsynons, data.frame(term, synons))

term <- "back pain"
synons <- c("back ache", "backpain", "lower back", "backache")
dfsynons <- rbind(dfsynons, data.frame(term, synons))

term <- "headache"
synons <- c("head ache", "head pain", "migraine")
dfsynons <- rbind(dfsynons, data.frame(term, synons))

term <- "uterus"
synons <- c("uterine")
dfsynons <- rbind(dfsynons, data.frame(term, synons))

term <- "temperature"
synons <- c("heat", "hot", "flush", "sweat")
dfsynons <- rbind(dfsynons, data.frame(term, synons))


term <- "fibroid"
synons <- c("fibrosis")
dfsynons <- rbind(dfsynons, data.frame(term, synons))

term <- "birth"
synons <- c("childbirth", "child", "kid", "partum", "baby", "children")
dfsynons <- rbind(dfsynons, data.frame(term, synons))

term <- "fall"
synons <- c("falling", "fell")
dfsynons <- rbind(dfsynons, data.frame(term, synons))


colnames(dfsynons) <- c("term", "synons")



## Process data:

dfQs <- as.data.frame(dfQs)

dfQs$allQs <- paste(dfQs$periodQs, dfQs$urinQs, dfQs$fibrQs,
                    dfQs$menopQs, dfQs$bowQs, dfQs$vagQs, dfQs$prolapseQs)



doQ <- function(df, model, Qcol, clinterms, fstword="pain") {
  
  clinterms <- unique(clinterms)
  
  df <- prep_Qcol(df, Qcol)
  
  df$wordMentioned <- ifelse(str_detect(df$allQs, fstword), "Mentioned", "Not Mentioned")

  # Pain
  dfPainM <- check_term(df, "pms", "wordMentioned", model)
  
  for(tm in unique(clinterms)){
    dfPainM <- rbind(dfPainM, check_term(df, tm, "wordMentioned", model))
  }
  
  # Race
  dfRace <- check_term(df, "pms", "Race", model)
  for(tm in clinterms){
  dfRace <- rbind(dfRace, check_term(df, tm, "Race", model))
  }
  
  #Age
  dfAge <- check_term(df, "pms", "Age", model)
  for(tm in clinterms){
    dfAge<- rbind(dfAge, check_term(df, tm, "Age", model))
  }
  
  #Job
  dfJob <- check_term(df, "pms", "Job_field", model)
  for(tm in clinterms){
    dfJob <- rbind(dfJob, check_term(df, tm, "Job_field", model))
  }
  
  #NC
  dfNC <- check_term(df, "pms", "Num_childbirths", model)
  for(tm in clinterms){
    dfNC <- rbind(dfNC, check_term(df, tm, "Num_childbirths", model))
  }
  
    
  gtRace <- test_term(dfRace,  clinterms, "Race", model, Qcol)
  gtAge <- test_term(dfAge,  clinterms, "Age", model, Qcol)
  gtJob <- test_term(dfJob,  clinterms, "Job_field", model, Qcol)
  gtNChild <- test_term(dfNC,  clinterms, "Num_childbirths", model, Qcol)
  gtNWordM <- test_term(dfPainM,  clinterms, "wordMentioned", model, Qcol)
  
  #tbres <- tbl_merge(list(gtRace, gtAge, gtJob, gtNChild))
  
  tbres <- list(gtRace, gtAge, gtJob, gtNChild, gtNWordM, df)
  
  return(tbres)
}


periodRes <- doQ(df=dfQs, 
                 model=model, 
                 Qcol="periodQs",
                 clinterms=clintermsPER
                 )

for(i in 1:5){
  gti <- periodRes[[i]]
  gt::gtsave(gti, filename = paste("results/period",i, ".pdf")
  )
}

urinRes <- doQ(df=dfQs, 
                 model=model, 
                 Qcol="urinQs",
                 clinterms=clintermsURIN)

for(i in 1:4){
  gti <- urinRes[[i]]
  gt::gtsave(gti, filename = paste("results/urine",i, ".pdf")
  )
}

fibrRes <- doQ(df=dfQs, 
               model=model, 
               Qcol="fibrQs",
               clinterms=clintermsFIBR)
for(i in 1:5){
  gti <- fibrRes[[i]]
  gt::gtsave(gti, filename = paste("results/fibr",i, ".pdf")
  )
}
vagRes <- doQ(df=dfQs, 
               model=model, 
               Qcol="vagQs",
               clinterms=clintermsVAG)
for(i in 1:5){
  gti <- vagRes[[i]]
  gt::gtsave(gti, filename = paste("results/vag",i, ".pdf")
  )
}

bowRes <- doQ(df=dfQs, 
              model=model, 
              Qcol="bowQs",
              clinterms=clintermsBOW)
for(i in 1:4){
  gti <- bowRes[[i]]
  gt::gtsave(gti, filename = paste("results/bow",i, ".pdf")
  )
}

prolRes <- doQ(df=dfQs, 
              model=model, 
              Qcol="prolapseQs",
              clinterms=clintermsPROL)
for(i in 1:4){
  gti <- prolRes[[i]]
  gt::gtsave(gti, filename = paste("results/prolapse",i, ".pdf")
  )
}


########### II. ADDITIONAL ANALYSES

###### 1. Word freq count overall

dfAll <- do.call("rbind", (list(periodRes[[5]],
                                urinRes[[5]],
                                bowRes[[5]],
                                vagRes[[5]],
                                prolRes[[5]],
                                fibrRes[[5]]
)))

tAll <- get_word_freq(dfAll, "allWfreq.txt")



###### 2. Word freq count by Qs

dfPer <- periodRes[[5]]
dfUri <- urinRes[[5]]
dfBow <- bowRes[[5]]
dfVag <- vagRes[[5]]
dfProl<- prolRes[[5]]
dfFibr <- fibrRes[[5]]



tPer <- get_word_freq(dfPer, "periodWfreq")
tUri <- get_word_freq(dfUri, "urineWfreq")
tBow <- get_word_freq(dfBow, "bowWfreq")
tVag <- get_word_freq(dfVag, "vagWfreq")
tProl <- get_word_freq(dfProl, "prolWfreq")
tFibr <- get_word_freq(dfFibr, "fibrWfreq")

###### 3. Word freq count by group

for(g in c("Race", "Num_childbirths", "Job_field","Age")) {
  
  t <- get_word_freq(dfPer, "periodWfreq", g)
  
  t <- get_word_freq(dfUri, "urineWfreq", g)
  t <- get_word_freq(dfBow, "bowWfreq", g)
  t <- get_word_freq(dfVag, "vagWfreq", g)
  t <- get_word_freq(dfProl, "prolWfreq", g)
  t <- get_word_freq(dfFibr, "fibrWfreq", g)
  
  t <- get_word_freq(dfAll, "allWfreq", g)
  
  # sent for all topics across groups as well (can consider also by topic but letc)
  ts1 <-syuzhet::get_nrc_sentiment(as.vector(t[[1]]$Var1))
  
  ts1$grp <- unique(df[g])[[1]][1]
  
  ts2 <-syuzhet::get_nrc_sentiment(as.vector(t[[2]]$Var1))
  ts2$grp <- unique(df[g])[[1]][2]
  
  ts <- rbind(ts1, ts2)
  
  gts <- gtsummary::tbl_summary(ts,
                         by = grp,
                         percent = "column", 
                         statistic = list(all_continuous() ~ "{median} ({p25}, {p75})",
                                          all_categorical() ~ "{n} ({p}%)"),
                         
  ) %>%
    gtsummary::modify_header(label = " ") %>%
    gtsummary::bold_labels() %>%
    gtsummary::modify_header(all_stat_cols() ~ "**{level}**, N = {n}") %>%
    gtsummary::add_p(test = all_categorical() ~ "fisher.test")
  
  gt::gtsave(as_gt(gts), filename = paste0("results/sent-by-",g, ".pdf"))
  
}

# wordclouds

pdf("results/wcAll.pdf")
  print(wordcloud::wordcloud(tAll$Var1, tAll$Freq, max.words = 250))
dev.off()

pdf("results/wcPer.pdf")
  print(wordcloud::wordcloud(tPer$Var1, tPer$Freq, max.words = 250))
dev.off()

pdf("results/wcVag.pdf")
  print(wordcloud::wordcloud(tVag$Var1, tVag$Freq, max.words = 250))
dev.off()

pdf("results/wcUri.pdf")
print(wordcloud::wordcloud(tUri$Var1, tUri$Freq, max.words = 250))
dev.off()

pdf("results/wcBow.pdf")
print(wordcloud::wordcloud(tBow$Var1, tBow$Freq, max.words = 250))
dev.off()

pdf("results/wcProl.pdf")
print(wordcloud::wordcloud(tProl$Var1, tProl$Freq, max.words = 250))
dev.off()



###### 4. Sentiment overall & Qs


s <- get_sents()

ggSent <- ggplot(s, aes(x=var, y=score, fill=var)) +
  geom_col() + 
  theme_bw() +
  theme(legend.title = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  facet_grid(~ sentiment)

pdf("results/sent.pdf", width=12, height=5)
  print(ggSent)
dev.off()

###### 5. Sentiment by group


##### 7. Conditional word analysis:

fstwords  <- c("pain", "cancer", "doctor", "bladder","birth", "maybe", "probably", "sometimes")

for( fstword in fstwords){
  fstword <- "bladder"
  allRes <- doQ(df=dfQs, 
                model=model, 
                Qcol="allQs",
                clinterms=c(clintermsPER, clintermsURIN, clintermsFIBR,
                            clintermsPROL, clintermsPROL, clintermsVAG,
                            clintermsBOW),
                fstword = fstword
  )
  
  gti <- allRes[[5]] %>%
    gt::tab_header(title = "Presence of terms, conditional on '",
                   str_to_title(fstword), 
                   "' being mentioned") %>%
    gt::gtsave(filename = paste("results/all_",
                                fstword , ".pdf")
    )
}

