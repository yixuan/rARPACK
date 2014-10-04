set.seed(123)
rands = runif(10000, -1, 1)
m = matrix(rands, 1000)
txt = apply(m, 1, paste, collapse = ", ")
txt = paste(txt, ",", sep = "")
txt[length(txt)] = gsub(",$", "", txt[length(txt)])
txt = c("const int rands_len = 10000;",
        "double rands[] = {",
        txt,
        "};")
write(txt, "rands.h")
