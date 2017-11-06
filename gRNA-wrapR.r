#gRNA-wrapR
#AUTOMATED PIPELINE FOR CRISPR/CAS9 gRNA AND PCR PRIMER DESIGN
#2017.3.9
#Gabriel Calvin

gRNA.wrapR = function(path, threshold, email, padding, cctop, mit.edu, entryDefault, entryMITSubmit, entryMITMonitor, entryMITGenebank, entryMITBlat, entryCCTOPSubmit, entryCCTOPMonitor, entryCCTOP.XLS, entryPostBlat, entryPrimer3, entryResults){
##############
#libraries to use for scraping
	if(require("rvest")){
		print("Library rvest is loaded correctly.")
	} else {
		print("Trying to install rvest.")
		install.packages("rvest")
		
		if(require(rvest)){
			print("Library rvest installed and loaded.")
		} else {
			stop("Could not install rvest.")
		}#end if-else
	}#end if-else
	
	if(require("RSelenium")){
		print("Library RSelenium is loaded correctly.")
	} else {
		print("Trying to install RSelenium.")
		install.packages("RSelenium")
		
		if(require(RSelenium)){
			print("Library RSelenium installed and loaded.")
		} else {
			stop("Could not install RSelenium.")
		}#end if-else
	}#end if-else
	
	require(rvest);	require(RSelenium)	
	
	if(missing(threshold)){
		threshold = 80
	}#end if
	
	check = FALSE
	while (check == FALSE) {
	
		if (anyNA(threshold) == TRUE) {
			threshold = ""
		}#end if
		
		#threshold
		#check that the threshold is an actual number/integer value between 0 and 100
		#sets to a default of 80%
		errorCheck = paste(" ", threshold, sep = "")
		if (nchar(errorCheck) == 1){
			threshold = 80
		} #end if
	
		if (threshold >= 0 & threshold <= 100){
			threshold = as.integer(threshold)
			print(paste("Threshold set to ", threshold, "%.", sep = ""))
			check = TRUE
		} else {
			threshold = 80
		}#end if-else
	}#end while
	#clear vars
	rm(errorCheck, check)
	
	#check that the padding is an actual number/integer value
	#sets to a default of 2500
	if(missing(padding)){
		padding = 2500
	}#end if
	
	if (anyNA(padding) == TRUE | is.numeric(padding) == FALSE){
		padding = 2500
	}#end if
	
	errorCheck = paste(" ", padding, sep = "")
	if (nchar(errorCheck) == 1){
		padding = 2500
	} #end if

	#padding
	print(paste("Padding set to ", padding, "bp.", sep = ""))	
	
	#clear vars
	rm(errorCheck)	
	
#swap backslash for forward slash and make it a path that R will recognize
	#path stays on memory
	path = gsub(pattern = "\\", replacement = "/", path, fixed = TRUE); #path	
	
#check that the path is an actual directory. bc trolls.
	isDir = strsplit(path, "/")[[1]]; #isDir
	isDir[1] = paste(isDir[1], "/", sep = ""); #isDir[1]
	
	n = 1
	while (n <= length(isDir)) {
		if (file.exists(isDir[1])) {
		
			if (n != 1) {
				isDir[n] = file.path(isDir[n-1], isDir[n])
				
				if (!file.exists(isDir[n])) {
					dir.create(file.path(isDir[n]))	
					print(paste("Created directory: ",isDir[n]))
				} #end if
			} #end if			
			
			#print(isDir[n])
			
			n = n + 1
		}#end if
		else {
			print("Directory not found.")
			n = length(isDir)+1
		}#end else
	}# end while-loop
	
	#clear variables
	rm(isDir, n)
	
#create directories for all the work results to be stored within
	#eight directories made so far
	#workFolders stays on memory
	workFolders = c("work", "geneSequences", "fastas", "compiledFastas", "MIT", "genebanks", "CCTOP", "zip", "XLS", "filteredGuides", "blat", "primer3", "results")
	if(!dir.exists(file.path(path, workFolders[1]))) {
		dir.create(file.path(path, workFolders[1]))
		dir.create(file.path(path, workFolders[1], workFolders[2]))
		dir.create(file.path(path, workFolders[1], workFolders[3]))
		dir.create(file.path(path, workFolders[1], workFolders[4]))
		dir.create(file.path(path, workFolders[1], workFolders[5]))
		dir.create(file.path(path, workFolders[1], workFolders[5], workFolders[6]))
		dir.create(file.path(path, workFolders[1], workFolders[7]))
		dir.create(file.path(path, workFolders[1], workFolders[7], workFolders[8]))
		dir.create(file.path(path, workFolders[1], workFolders[7], workFolders[9]))
		dir.create(file.path(path, workFolders[1], workFolders[10]))
		dir.create(file.path(path, workFolders[1], workFolders[11]))
		dir.create(file.path(path, workFolders[1], workFolders[12]))	
	} else {
		print("A work folder already exists for this directory.")
	}#end if-else
	
	if(!dir.exists(file.path(path,workFolders[13]))){
		dir.create(file.path(path,workFolders[13]))
	} else {
		print("A results folder already exists for this directory.")	
	}#end if-else
##############

##########
#switching off modules to allow different entry points in the pipeline
if(missing(entryDefault)){
	entryDefault = TRUE
}#end if

	#MIT entry points
	if(missing(entryMITSubmit)){
		entryMITSubmit = FALSE
	}#end if
	if(missing(entryMITMonitor)){
		entryMITMonitor = FALSE
	}#end if
	if(missing(entryMITGenebank)){
		entryMITGenebank = FALSE
	}#end if
	if(missing(entryMITBlat)){
		entryMITBlat = FALSE
	}#end if

	#CCTOP entry points
	if(missing(entryCCTOPSubmit)){
		entryCCTOPSubmit = FALSE
	}#end if
	if(missing(entryCCTOPMonitor)){
		entryCCTOPMonitor = FALSE
	}#end if
	if(missing(entryCCTOP.XLS)){
		entryCCTOP.XLS = FALSE
	}#end if

	#BACKEND entry points
	if(missing(entryPostBlat)){
		entryPostBlat = FALSE
	}#end if
	if(missing(entryPrimer3)){
		entryPrimer3 = FALSE
	}#end if
	if(missing(entryResults)){
		entryResults = FALSE
	}#end if
	
	if(entryMITSubmit == TRUE | entryMITMonitor == TRUE | entryMITGenebank == TRUE | entryMITBlat == TRUE | entryCCTOPSubmit == TRUE | entryCCTOPMonitor == TRUE | entryCCTOP.XLS == TRUE | entryPostBlat == TRUE | entryPrimer3 == TRUE | entryResults == TRUE){
		entryDefault = FALSE
	}#end if
############

#############
#decide between cctop or crispr.mit.edu
	check = FALSE
	if(missing(cctop)){
		cctop = TRUE
	}#end if
	
	if(missing(mit.edu)){
		mit.edu = FALSE
	}#end if
	
	while(check == FALSE){
		if (cctop == FALSE & mit.edu == TRUE){
			print("crispr.mit selected")
			check = TRUE
		} else if (cctop == TRUE & mit.edu == FALSE){
			print("CCTOP selected")
			check = TRUE
		} else if (cctop == "" & mit.edu == "") {
			cctop = TRUE
			mit.edu = FALSE				
		} else if (cctop == TRUE & mit.edu == TRUE) {
			stop("Error: Only one prediction software can be selected")
		}#end if-else
	}#end while
	
	rm(check); gc()	
#############

#############
#save options as a file
	dateTime = strsplit(as.character(Sys.time()), split = " ")[[1]]
	dateTime = paste(dateTime[1], dateTime[2], sep = "_"); #dateTime
	dateTime = gsub(dateTime, pattern = "-|:", replacement = ""); #dateTime
	write(c(path, threshold, email, padding, cctop, mit.edu, entryDefault, entryMITSubmit, entryMITMonitor, entryMITGenebank, entryMITBlat, entryCCTOPSubmit, entryCCTOPMonitor, entryCCTOP.XLS, entryPostBlat, entryPrimer3, entryResults), file.path(path, paste("selectedOptions_", dateTime, ".txt", sep = "")))
	rm(dateTime); gc()
	print(paste("Job started at:  ", Sys.time(), sep = ""))
#############

##############
#Scan a file from the working folder. This first file should contain the names of the genes to get information from Gene Cards, genecards.org as a .txt file
if(entryDefault == TRUE){
#it takes the website and combines it with the name of the gene.
#it then scrapes the combined url for the information of interest:
	#chromosome
	#start
	#end
	#orientation
	#size in bp
#this information is combined in the form of chromosome:start-end
	#whole gene coordinates are then used in another section that will download the sequence from the ucsc browser.
#coordinates for upstream and downstream areas of interest may have to be entered into a file manually, which will also contain the information for the gene.

#take values from a file and put them together with the genecards website to generate urls to gene pages, which will be scraped
	
	#scan the initial document
	if(file.exists(file.path(path, "genes0.txt"))){
		genes0 = list.files(path, full.names = TRUE, pattern = "genes0.txt")
		#geneNameFromFile stays on memory
		geneNameFromFile = scan(genes0, what = " ", quiet = TRUE); #geneNameFromFile
	}else {
		stop("No initial genes0.txt file found")
	}#end if
	
	#clear variables
	rm(genes0)
###########	

######
#sift the entries in the genes file, and detect whether it is either 1) a gene name; 2) a sequence position; 3)a refseq; 4) something other.

	#detect based on ":" for position
		if (length(grep(geneNameFromFile, pattern = ":")) >= 1){
			#sift into position variable
			chrCoordinates = geneNameFromFile[grep(geneNameFromFile, pattern = ":")]
			#remove from geneNameFromFile
			geneNameFromFile = geneNameFromFile[-c(grep(geneNameFromFile, pattern = ":"))]
			
			for(a0 in 1:length(chrCoordinates)){
				chrCoordinates.chr = strsplit(chrCoordinates[a0], split = ":")[[1]]; #chrCoordinates.chr
				chrCoordinates.start = strsplit(chrCoordinates.chr[2], split = "-")[[1]][1]; #chrCoordinates.start
				chrCoordinates.end = strsplit(chrCoordinates.chr[2], split = "-")[[1]][2]; #chrCoordinates.end
				
				chrCoordinates.chr = chrCoordinates.chr[1]; #chrCoordinates.chr
				
				if(a0 == 1){
					chrCoordinates.matrix = matrix(ncol = 3, data = c(chrCoordinates.chr, chrCoordinates.start, chrCoordinates.end), byrow = TRUE)
				} else {
					chrCoordinates.matrix = rbind(chrCoordinates.matrix, chrCoordinates.chr, chrCoordinates.start, chrCoordinates.end)
				}#end if-else
			}#end a0 for
			
			write.csv(chrCoordinates.matrix, file.path(path, "chrCoordinates.csv"), row.names = FALSE, quote = FALSE)
		#clear vars
		rm(chrCoordinates, chrCoordinates.chr, chrCoordinates.start, chrCoordinates.end, chrCoordinates.matrix, a0)
		} #end if
		
		
		
	#detect based on "_" for refseqs
		if (length(grep(geneNameFromFile, pattern = "_")) >= 1){
			#check against refseq prefixes to determine if genomic or nongenomic
			refSeqPrefixes.genomic = c("AC_", "NC_", "NG_", "NT_", "NW_", "NZ_")
			refSeqPrefixes.nongenomic = c("NM_", "NR_", "XM_", "XR_", "AP_", "NP_", "YP_", "XP_", "WP_")
			
			#sift into refseq variable
			refSeqID = geneNameFromFile[grep(geneNameFromFile, pattern = "_")]
			#remove from geneNameFromFile	
			geneNameFromFile = geneNameFromFile[-c(grep(geneNameFromFile, pattern = "_"))]
			
			#dummy vector
			refSeqFilter.nongenomic = ""
			
			for (a1 in 1:length(refSeqID)){
				if (length(which(substr(refSeqID[a1], start = 1, stop = 3) == refSeqPrefixes.genomic)) >= 1){
					if (a1 == 1){
						refSeqFilter = refSeqID[a1]
					} else {
						refSeqFilter[a1] = refSeqID[a1]
					} #end if-else
				} else {
					refSeqFilter.nongenomic = c(refSeqFilter.nongenomic, refSeqID[a1])
				}#end if-else			
			} #end a1 for
			#refSeqFilter
			
			if (length(refSeqFilter.nongenomic) >= 1) {
				write(refSeqFilter.nongenomic, file.path(path, "rejectedEntries.txt"))
			} #end if
			
			rm(refSeqPrefixes.genomic, refSeqPrefixes.nongenomic, refSeqFilter.nongenomic, refSeqID, a1)

		} #end if
		
	#set aside anything that isn't one of those and say that it couldn't figure it out as a report file

	#make geneNameFromFile1 variable with new list which should contain only gene names now
	geneNameFromFile1 = geneNameFromFile

######

##############	
	#RSelenium
	#open the browser and navigate to genecards.org
	RSDRIVER = rsDriver(port = 4444L, browser = "chrome", check = TRUE, verbose = FALSE)
	RS.browser = RSDRIVER$client

	for(b0 in 1:length(geneNameFromFile)){
		geneNameFromFile[b0] = paste("http://www.genecards.org/cgi-bin/carddisp.pl?gene=", geneNameFromFile[b0], "&keywords=", geneNameFromFile[b0], sep = "")
		#print(geneNameFromFile[b0])
		
		RS.browser$navigate(geneNameFromFile[b0])
		Sys.sleep(5)
		#genecards.org information
		#css selector for the chromosome number
		# "css selector', '#genomic_location > div:nth-child(3) > div:nth-child(3) > div:nth-child(2) > dl:nth-child(1) > dd:nth-child(2)"
		gene.chrNum = RS.browser$findElement(using = 'css selector', '#genomic_location > div:nth-child(3) > div:nth-child(3) > div:nth-child(2) > dl:nth-child(1) > dd:nth-child(2)'); #gene.chrNum
		gene.chrNum = gene.chrNum$getElementText()[[1]]; #gene.chrNum
		gene.chrNum = paste("chr", gene.chrNum, sep = ""); #gene.chrNum		
		
		#css selector for start coordinates
		# "css selector', '#genomic_location > div:nth-child(3) > div:nth-child(3) > div:nth-child(2) > dl:nth-child(2) > dd:nth-child(2)"
		gene.start = RS.browser$findElement(using = 'css selector', '#genomic_location > div:nth-child(3) > div:nth-child(3) > div:nth-child(2) > dl:nth-child(2) > dd:nth-child(2)'); #gene.start
		gene.start = gene.start$getElementText()[[1]]; #gene.start
		gene.start = strsplit(gene.start, split = " ")[[1]][1]; #gene.start
		gene.start = gsub(x = gene.start, pattern = ",", replacement = ""); #gene.start		
			
		#css selector for end coordinates		
		# "#genomic_location > div:nth-child(3) > div:nth-child(3) > div:nth-child(2) > dl:nth-child(2) > dd:nth-child(4)"		
		gene.end = RS.browser$findElement(using = 'css selector', '#genomic_location > div:nth-child(3) > div:nth-child(3) > div:nth-child(2) > dl:nth-child(2) > dd:nth-child(4)'); #gene.end
		gene.end = gene.end$getElementText()[[1]]; #gene.end
		gene.end = strsplit(gene.end, split = " ")[[1]][1]; #gene.end
		gene.end = gsub(x = gene.end, pattern = ",", replacement = ""); #gene.end
				
		#css selector for size
		# "dl.dl-inline:nth-child(3) > dd:nth-child(2)"		
		gene.size = RS.browser$findElement(using = 'css selector', 'dl.dl-inline:nth-child(3) > dd:nth-child(2)'); #gene.size
		gene.size = gene.size$getElementText()[[1]]; #gene.size
		gene.size = strsplit(gene.size, " ")[[1]][1]; #gene.size
		gene.size = gsub(x = gene.size, pattern = ",", replacement = ""); #gene.size
		
		#css selector for orientation
		# "dl.dl-inline:nth-child(3) > dd:nth-child(4)"		
		gene.orientation = RS.browser$findElement(using = 'css selector', 'dl.dl-inline:nth-child(3) > dd:nth-child(4)'); #gene.orientation
		gene.orientation = gene.orientation$getElementText()[[1]]; #gene.orientation
		gene.orientation = strsplit(gene.orientation, " ")[[1]][1]; #gene.orientation
	
		#css selector for the link to the genecards ucsc browser custom track
		# "#genomic_location > div:nth-child(3) > div:nth-child(4) > div:nth-child(2) > div:nth-child(1) > a:nth-child(1)"		
		gene.customTrackLink = RS.browser$findElement(using = 'css selector', '#genomic_location > div:nth-child(3) > div:nth-child(4) > div:nth-child(2) > div:nth-child(1) > a:nth-child(1)'); #gene.customTrackLink
		gene.customTrackLink = gene.customTrackLink$getElementAttribute(attr = 'href')[[1]]; #gene.customTrackLink
		
		#check that the link for the ucsc browser is the correct position
		coordCheck = strsplit(gene.customTrackLink, ":")[[1]]; #coordCheck
		coordCheck = coordCheck[length(coordCheck)]; #coordCheck
		coordCheck = strsplit(coordCheck, "-")[[1]]; #coordCheck
		
		if (coordCheck[1] != gene.start){
			gene.customTrackLink = gsub(gene.customTrackLink, pattern = coordCheck[1], replacement = gene.start); #gene.customTrackLink
		} #end if
		
		if (coordCheck[2] != gene.end){
			gene.customTrackLink = gsub(gene.customTrackLink, pattern = coordCheck[2], replacement = gene.end); #gene.customTrackLink
		} #end if
		
		if(geneNameFromFile[b0] == geneNameFromFile[1]) {
			geneName.matrix = matrix(ncol = 8, data = c(geneNameFromFile1[b0], geneNameFromFile[b0], gene.chrNum, gene.start, gene.end, gene.size, gene.orientation, gene.customTrackLink)); #geneName.matrix
			
			} else {
			geneName.matrix = rbind(geneName.matrix, c(geneNameFromFile1[b0], geneNameFromFile[b0], gene.chrNum, gene.start, gene.end, gene.size, gene.orientation, gene.customTrackLink)); #geneName.matrix
		} #end if-else
	} #end b0 for-loop
	
#save the genNameFromFile in the work workFolders[1]
	write(geneNameFromFile, file.path(path, workFolders[1], "genecards_geneURLs.txt"))
	
#save the geneName.matrix in the work workFolders[1]
	write.csv(geneName.matrix, file.path(path, workFolders[1], "genecards_geneInformation.csv"), row.names = FALSE, quote = FALSE)

#close browser and quit server
	RS.browser$close()
	RSDRIVER$server$stop()
	
#remove variables for good housekeeping
	rm(gene.chrNum, gene.start, gene.end, gene.size, gene.orientation, gene.customTrackLink, coordCheck, b0)
	#geneName.matrix stays on memory
	
#############

##########
############IN DEVELOPMENT###########
#obtain genomic coordinates from refseqID
#for (g0 in 1:length(refSeqFilter)){
	#navigate browser to refseq website
	#RS.browser$navigate("https://www.ncbi.nlm.nih.gov/nuccore")

	#css selector for search bar
	# "#term"
	#refSeq.search = RS.browser$findElement(using = 'css', "#term")
	#refSeq.search$clickElement()
	#refSeq.search = RS.browser$sendKeysToActiveElement(list(refSeqFilter[g0]))
	#refSeq.search = RS.browser$sendKeysToActiveElement(list("NC_000023.11"))
	
	#css selector for submit
	# "#search"
	#refSeq.submit = RS.browser$findElement(using = 'css', "#search")
	#refSeq.submit$clickElement()
	
	#css selector for fasta
	# "#ReportShortCut6"
	#refSeq.fasta = RS.browser$findElement(using = 'css', "#ReportShortCut6")
	#refSeq.fasta$clickElement()
	#refSeq.fasta = RS.browser$getPageSource()[[1]]
	
	#} #end for g0

############

#############
#use the coordinates contained in the genecards_geneInformation.csv file to get sequence from the UCSC browser as fasta files.
#the fasta files will be contained in the working folder, and will be chopped into smaller pieces of sequence, and then saved as part of a larger fasta file.
	#~28 sequences fit under the size limit for crispr.mit.edu
	
#picks up the 8th column in the file, which corresponds to the link scraped from genecards.org
	for (c0 in 1:length(geneNameFromFile1)) {
		ucsc.geneBrowserURL = read_html(geneName.matrix[c0, 8])
				
		ucsc.viewDNA = html_attr(html_nodes(ucsc.geneBrowserURL, css = "#dnaLink"), "href"); ucsc.viewDNA
		ucsc.viewDNA = substr(ucsc.viewDNA, 4, nchar(ucsc.viewDNA)); ucsc.viewDNA
		ucsc.viewDNA = paste("genome.ucsc.edu/", ucsc.viewDNA, sep = ""); ucsc.viewDNA
		
		ucsc.viewDNA.session = html_session(ucsc.viewDNA)
		
		ucsc.viewDNA.form = html_form(ucsc.viewDNA.session)[[1]]
		ucsc.viewDNA.form$fields = ucsc.viewDNA.form$fields[-c(15,19)]
		ucsc.viewDNA.form = set_values(form = ucsc.viewDNA.form, hgSeq.padding5 = padding, hgSeq.padding3 = padding)
		ucsc.viewDNA.sequenceURL = submit_form(ucsc.viewDNA.session, ucsc.viewDNA.form); #ucsc.viewDNA.sequenceURL
		
		#file is downloaded to the geneSequences workFolders[2]
		download.file(ucsc.viewDNA.sequenceURL$url, file.path(path, workFolders[1], workFolders[2], paste("UCSC_", geneNameFromFile1[c0], "_sequence_padded", padding,"bp_wholeGene.fa", sep = "")))
	}#end c0 for-loop
	
	#clear variables
	rm(c0, ucsc.geneBrowserURL, ucsc.viewDNA, ucsc.viewDNA.session, ucsc.viewDNA.form, ucsc.viewDNA.sequenceURL)

############	
#uses the chromosome position to obtain the sequence
	if(file.exists(file.path(path, "chrCoordinates.csv"))){
		chrCoordinates.filePath = list.files(path, pattern = "chrCoordinates.csv", full.names = TRUE)
		
		chrCoordinates.matrix = read.csv(chrCoordinates.filePath, header = TRUE); chrCoordinates.matrix
		chrCoordinates = paste(paste(chrCoordinates.matrix[,1], chrCoordinates.matrix[,2], sep = ":"), chrCoordinates.matrix[,3], sep = "-")
		RSDRIVER = rsDriver(port = 4444L, browser = "chrome", check = TRUE, verbose = FALSE)
		RS.browser = RSDRIVER$client

		for(c5 in 1:length(chrCoordinates)){	
			RS.browser$navigate("http://genome.ucsc.edu/cgi-bin/hgGateway")
			
			Sys.sleep(1)

			#css selector for position input for ucsc
			#"#positionInput"
			ucsc.input = RS.browser$findElement(using = 'css', "#positionInput")
			ucsc.input$clickElement()
			ucsc.input = RS.browser$sendKeysToActiveElement(list(chrCoordinates[c5]))
			
			#css selector for submit button
			#".jwGoButton"
			ucsc.submit = RS.browser$findElement(using = 'css', '.jwGoButton')
			ucsc.submit$clickElement()
			
			#css selector for viewDNA
			#"#dnaLink"
			ucsc.viewDNA = RS.browser$findElement(using = 'css', "#dnaLink")
			RS.browser$navigate(ucsc.viewDNA$getElementAttribute(attr = 'href')[[1]])
			
			#xpath selector for 5' padding in getDNA
			# '//*[@id="hgSeq.padding5"]'
			ucsc.getDNA5 = RS.browser$findElement(using = 'xpath', '//*[@id="hgSeq.padding5"]')
			ucsc.getDNA5$setElementAttribute(attributeName = 'value', padding)
			
			#xpath for 3' padding in getDNA
			#'//*[@id="hgSeq.padding3"]'
			ucsc.getDNA3 = RS.browser$findElement(using = 'xpath', '//*[@id="hgSeq.padding3"]')
			ucsc.getDNA3$setElementAttribute(attributeName = 'value', padding)
			
			#css for submit
			#"#firstSection > table:nth-child(1) > tbody:nth-child(1) > tr:nth-child(1) > td:nth-child(1) > table:nth-child(1) > tbody:nth-child(1) > tr:nth-child(1) > td:nth-child(1) > table:nth-child(2) > tbody:nth-child(1) > tr:nth-child(2) > td:nth-child(2) > form:nth-child(5) > p:nth-child(30) > input:nth-child(1)"
			ucsc.submit = RS.browser$findElement(using = 'css', "#firstSection > table:nth-child(1) > tbody:nth-child(1) > tr:nth-child(1) > td:nth-child(1) > table:nth-child(1) > tbody:nth-child(1) > tr:nth-child(1) > td:nth-child(1) > table:nth-child(2) > tbody:nth-child(1) > tr:nth-child(2) > td:nth-child(2) > form:nth-child(5) > p:nth-child(30) > input:nth-child(1)")
			ucsc.submit$clickElement()
			
			Sys.sleep(10)
			
			ucsc.viewDNA.sequenceURL = RS.browser$getCurrentUrl()[[1]]
			#file is downloaded to the geneSequences workFolders[2]
			download.file(ucsc.viewDNA.sequenceURL, file.path(path, workFolders[1], workFolders[2], paste("UCSC_", gsub(chrCoordinates[c5], pattern = ":", replacement = "-"), "_sequence_padded", padding,"bp_wholeGene.fa", sep = "")))
			
			Sys.sleep(1)
		}#end c5 for
		
		RS.browser$close()
		RSDRIVER$server$stop()
	}#end if
############
	
	#ucsc.viewDNA.filePath stays in memory for a bit
	#picks up whole gene sequences in the geneSequences workFolders[2]
	ucsc.viewDNA.filePath =  list.files(file.path(path, workFolders[1], workFolders[2]), full.names = TRUE, pattern = "wholeGene.fa")
	
	for (c1 in 1:length(ucsc.viewDNA.filePath)){
		ucsc.viewDNA.file = scan(ucsc.viewDNA.filePath[c1], what = "", sep = "\n", quiet = TRUE)
		
		check = FALSE
		
		while (check == FALSE) {
			if(substr(ucsc.viewDNA.file[1], 1, 5) == ">hg38") {
				check = TRUE
			}#end if	
			else {
				ucsc.viewDNA.file = ucsc.viewDNA.file[-1]
			}#end else
		}#end while-loop
		
		ucsc.viewDNA.fileName = strsplit(ucsc.viewDNA.filePath[c1], "/")[[1]]
		ucsc.viewDNA.fileName = strsplit(ucsc.viewDNA.fileName[length(ucsc.viewDNA.fileName)], "_")[[1]]
		ucsc.viewDNA.file[1] = paste(">", ucsc.viewDNA.fileName[2], sep = "")			
		
		ucsc.viewDNA.file = ucsc.viewDNA.file[-c(which(ucsc.viewDNA.file == "</PRE>"):length(ucsc.viewDNA.file))]
		write(ucsc.viewDNA.file, ucsc.viewDNA.filePath[c1])
	} #end c1 for-loop
	
	#clear vars
	rm(c1, ucsc.viewDNA.fileName, ucsc.viewDNA.file)
	
	for (c2 in 1:length(ucsc.viewDNA.filePath)){
		#while-loop for splitting fastas
		fileToSplit = scan(ucsc.viewDNA.filePath[c2], what = "", sep = "\n", quiet = TRUE)
		c2.i = 2
		c2.1 = 2
		c2.2 = 6
		c2.Ticker = 1
		
		check = FALSE
		fastaSeqName = fileToSplit[1]; fastaSeqName = substr(fastaSeqName, 2, nchar(fastaSeqName))
	
		while (c2.i <= length(fileToSplit)){
			fileToSplit.part = fileToSplit[c2.1:c2.2]; #print(fileToSplit.part)
			
			#file saved to the fastas workFolders[3]
			write(c(paste(">", fastaSeqName, "_", c2.Ticker, sep = ""), fileToSplit.part), file.path(path, workFolders[1], workFolders[3], paste(fastaSeqName, "_", c2.Ticker, ".fa", sep = "")))
			
			#start next sequence at last line of the previous to achieve tiling
			c2.1 = c2.2
			c2.2 = c2.2 + 4
			c2.Ticker = c2.Ticker + 1
			c2.i = c2.2
			
			if (check == FALSE){
				if(c2.2 > length(fileToSplit)){
					c2.2 = length(fileToSplit)
					c2.i = c2.2
					check = TRUE
				}#end if
			}#end if
		}#end while-loop
	}#end c2 for-loop
	
	#clear vars
	rm(c2, fileToSplit, c2.i, c2.1, c2.2, c2.Ticker, fastaSeqName, fileToSplit.part, check, ucsc.viewDNA.filePath)
		
	#compiles fastas into larger fasta files, of 28 sequences each, which will fit into crispr.mit.edu, since their size limit is 11kb
	#picks up fasta files in the fastas workFolders[3]
	files.fa = list.files(file.path(path, workFolders[1], workFolders[3]), full.names = TRUE, pattern = ".fa"); #files.fa

	c4.i = 1
	c4.1 = 1
	c4.2 = 28
	check = FALSE
	
	while(c4.i <= ((length(files.fa)%/%28)+1)){
		compiledFa = scan(files.fa[c4.1], what = "", quiet = TRUE); #compiledFa
		lastUnderscore = which(strsplit(compiledFa[1], "")[[1]] == "_")
		seqName = substr(compiledFa[1], 2, lastUnderscore[length(lastUnderscore)]-1); #print(seqName)	
		for (c4.0 in (c4.1 + 1):c4.2){
			fa.filesToScan = scan(files.fa[c4.0], what = "", quiet = TRUE)
			compiledFa = c(compiledFa, fa.filesToScan)		
		}#end for-loop
	
		c4.1 = c4.1 + 28
		c4.2 = c4.2 + 28
		
		#saves files to compiledFastas workFolders[4]
		write(compiledFa, file.path(path, workFolders[1], workFolders[4], paste(seqName, c4.i, "compiled.fa", sep = "_")))
		
		c4.i = c4.i + 1
		
		if (check == FALSE){
			if(c4.2 > length(files.fa)){
				c4.2 = length(files.fa)
				check = TRUE
			}#end if
		}#end if		
	}#end while-loop
	
	#clear vars	
	rm(files.fa, c4.i, c4.1, c4.2, c4.0, compiledFa, seqName, fa.filesToScan, lastUnderscore, entryDefault)
	gc()
############

entryMITSubmit = TRUE
entryCCTOPSubmit = TRUE
}#end entryDefault if
#############

#############
#MIT
	#take the fasta files that were generated, which contain smaller bits of the sequence from the original fastas, and submit them as a job to the crispr.mit.edu batch page.
	#once a batch is submitted, url for the job status page to be stored in another file.
if(cctop == FALSE & mit.edu == TRUE) {
	if(entryMITSubmit == TRUE){
	#list the compiledFa.fa files in the compiledFastas workFolders[4]
	files.compiledFa = list.files(file.path(path, workFolders[1], workFolders[4]), full.names = TRUE, pattern = "compiled.fa")
	
	if(length(files.compiledFA) >= 1){
	#RSelenium
	#open the browser and navigate to crispr.mit.edu
		RSDRIVER = rsDriver(port = 4444L, browser = "chrome", check = TRUE, verbose = FALSE)
		RS.browser = RSDRIVER$client

	#d0 for-loop 
	for (d0 in 1:length(files.compiledFa)){
		RS.browser$navigate("http://crispr.mit.edu")
		Sys.sleep(10)
		
		#click over to the batch entry page
		mit.batchEntry = RS.browser$findElement(using = 'css selector', ".nav-tabs > li:nth-child(2)")
		mit.batchEntry$clickElement()
		
		#give the select-file element the fasta
		mit.fastaSelect = RS.browser$findElement(using = 'css selector', "#file-input-2")
			tempFasta = tempfile(fileext = ".fa")
			fastaFile = scan(file = files.compiledFa[d0], what = "", sep = "\n", quiet = TRUE)
			write(fastaFile, tempFasta)
		mit.fastaSelect$sendKeysToElement(list(tempFasta))
		
		#give the email address
		mit.emailEnter = RS.browser$findElement(using = 'css selector', "#email-address-2")
		mit.emailEnter$sendKeysToElement(list(email))
		
		#submit the form
		mit.clickSubmit = RS.browser$findElement(using = 'css selector', "div.submit-group:nth-child(4) > div:nth-child(4) > input:nth-child(1)")
		mit.clickSubmit$clickElement()
		
		#capture the url of the job status page
		Sys.sleep(30) #needs to wait an amount of time for the job to refresh. maybe wait loop?
		if (d0 == 1){
			mit.url = RS.browser$getCurrentUrl()[[1]]; #mit.url
		} else {
			mit.url = c(mit.url, RS.browser$getCurrentUrl()[[1]])
		} #end if-else
	} #end d0 for-loop
	
	#close the browser and stop the server
	RS.browser$close()
	RSDRIVER$server$stop()
	
	#save jobStatusURLs to the MIT workFolders[5]
	write(mit.url, file.path(path, workFolders[1], workFolders[5], "mit_jobStatusURLs.txt"))

	#clear vars
	rm(files.compiledFa, mit.batchEntry, mit.fastaSelect, tempFasta, fastaFile, mit.emailEnter, mit.clickSubmit, mit.url, d0, entryMITSubmit)
	gc()
	
	entryMITMonitor = TRUE
	
	} else {
		stop("No compiled.fa files were found")
	}#end if
	
	}#end entryMITSubmit if
############

#############
#scan a file that contains all the urls for the job status webpages, and then loops for each one.	
#scrape the job status webpage from crispr.mit.edu
if(entryMITMonitor == TRUE){
	gb.geneList = scan(file = file.path(path, workFolders[1], workFolders[5], "mit_jobStatusURLs.txt"), what = "", sep = "\n", quiet = TRUE); #gb.geneList
	
	if(length(gb.geneList) >= 1){
	for (d1 in 1:length(gb.geneList)) {
		gb.gene = scan(file = gb.geneList[d1], what = "", quiet = TRUE)
		
		#parse until "sequence"
		check = FALSE; #check
		 while (check == FALSE) {
			if (isTRUE(gb.gene[1] == "sequence")) {
				check = TRUE; #check
			} else {
				gb.gene = gb.gene[-1]; #gb.gene
				} #end if-else		
		} #end while-loop

		#make a vector with the urls for the jobs
		gb.gene = gb.gene[-c(which(gb.gene == ":"))]; gb.gene = gb.gene[-c(which(gb.gene == ","))]
		gb.gene = gb.gene[-c(which(gb.gene == "original_filename"):length(gb.gene))]; #gb.gene
		gb.gene = gb.gene[-c((length(gb.gene)-9):length(gb.gene))]; #gb.gene

		if (length(gb.gene)%%63 == 0) {
			gb.geneMatrix = matrix(ncol = 63, nrow = length(gb.gene)%/%63, data = gb.gene, byrow = T)	
		} #end if

		#create the URL for the genebank, this is column 45 in the gb.geneMatrix
		gb.geneMatrix[,45] = paste("http://crispr.mit.edu/export/guides_gb/", gb.geneMatrix[,45], sep = ""); #gb.geneMatrix[,45]
		#get the names of the files used for the genebank URL, this is column 53 in the gb.geneMatrix
		
		if (d1 == 1){
			genes.url = gb.geneMatrix[,45]
			genes.names = gb.geneMatrix[,53]
		} else {
			genes.url = c(genes.url, gb.geneMatrix[,45])
			genes.names = c(genes.names, gb.geneMatrix[,53])
		}#end if-else
	}#end d1 for-loop
	
	#make a vector containing the job URL
	genes.jobURL = gsub(genes.url, pattern = "export/guides_gb", replacement = "job"); #genes.jobURL
	
	#save the guide URLs to MIT workFolders[5]	
	gb.geneMatrix = matrix(ncol = 3, data = c(genes.names, genes.jobURL, genes.url), byrow = FALSE)
	write.csv(gb.geneMatrix, file.path(path, workFolders[1], workFolders[5], paste("mit_genebankURLs.csv", sep = "")), row.names = FALSE)	
	
	#make a dummy vector
	genes.filepaths = 1:length(genes.names); #genes.filepaths
		
	#download the genebank files into the geneBanks workFolders[6]
	#monitor for readiness
	for (d2 in 1:length(genes.names)){
		check = FALSE
		while (check == FALSE){
		checkExists = scan(genes.jobURL[d2], what = "", quiet = TRUE, sep = ""); #checkExists
			#checks for the "files_ready" value to be true or false
			if(checkExists[76] != "false,"){
				check = TRUE
			} else {
				print(paste(Sys.time(), "  :  Nothing just yet...", sep = ""))
				Sys.sleep(180)
			}#end if-else
		}#end while
		
		genes.filepaths[d2] = file.path(path, workFolders[1], workFolders[5], workFolders[6], paste(genes.names[d2], ".gb", sep = ""))
		download.file(genes.url[d2], genes.filepaths[d2])	
	}#end d2 for-loop
	
	rm(entryMITMonitor)
	gc()
	
	entryMITGenebank = TRUE
	} else {
		stop("No jobStatusURLs.txt file was found")
	}#end if
}#end entryMITMonitor if

if(entryMITGenebank == TRUE){
	#make a list of the filepaths in the geneBanks workFolders[6]
	gb.filepaths = list.files(file.path(path, workFolders[1], workFolders[5], workFolders[6]), full.names = TRUE, pattern = ".gb"); #gb.filepaths
	
	if(length(gb.filepaths) >= 1){
	#dummy vectors
	name.previous = ""; #name.previous
	compiled.matrix = matrix(nrow = 1, ncol = 5); #compiled.matrix

	for (d3 in 1:length(gb.filepaths)){
		genebank = scan(file = gb.filepaths[d3], what = " ", quiet = TRUE); #genebank
		name = gsub(pattern = ",", replacement = "", x = genebank[13]); #name
		lastUnderscore = which(strsplit(name, "")[[1]] == "_")
		name = substr(name, 1, lastUnderscore[length(lastUnderscore)]-1); #name	
		
		#parse until "protein_bind"
		check = FALSE; #check

		while (check == FALSE) {
			if (isTRUE(genebank[1] == "protein_bind")) {
				check = TRUE
				} else {
				genebank = genebank[-1]
				}#end if-else			
		} #end while-loop

		#genebank
		genebank = genebank[-c(which(genebank == "ORIGIN"):length(genebank))]; #genebank

		genebank.matrix = matrix(ncol = 20, nrow = length(genebank)%/%20, data = genebank, byrow = T); #genebank.matrix
		genebank.matrix[,1] = name
		genebank.matrix[,2] = as.integer(gsub(pattern = "%", replacement = "", x = genebank.matrix[,15])); #genebank.matrix
		genebank.matrix[,4] = paste(genebank.matrix[,4], genebank.matrix[,5], genebank.matrix[,6]); #genebank.matrix
		genebank.matrix[,8] = paste(genebank.matrix[,8], genebank.matrix[,9], genebank.matrix[,10], genebank.matrix[,11], genebank.matrix[,12]); #genebank.matrix
		genebank.matrix = genebank.matrix[,-c(3,5,6,7,9:18,20)]; #genebank.matrix
		######THIS IS WHERE GUIDES ARE FILTERED ACCORDING TO THE THRESHOLD SET AT THE BEGINNING######		
		genebank.matrix = genebank.matrix[-c(which(genebank.matrix[,2] < threshold)),]; #genebank.matrix
	
	if(length(genebank.matrix) != 0){
		############
		#reverse complement
			guideList = genebank.matrix[,5]; guideList
			
			#remove PAM NGG
			guideList = substr(guideList, 1, (nchar(guideList)-3)); guideList
			
			for(d3.1 in 1:length(guideList)){
				#add 'G' for U6 promotor if necessary 
				if (substr(guideList[d3.1], 1, 1) != "G"){
					guideList[d3.1] = paste("G", guideList[d3.1], sep = "")
				} #end if
				
				guideList.split = strsplit(guideList[d3.1], "")[[1]]; guideList.split

				for (d3.1.1 in 1:length(guideList.split)){
					if(guideList.split[d3.1.1] == "G"){
						guideList.split[d3.1.1] = "C"
					} else if(guideList.split[d3.1.1] == "C") {
						guideList.split[d3.1.1] = "G"
					} else if(guideList.split[d3.1.1] == "A") {
						guideList.split[d3.1.1] = "T"
					} else if(guideList.split[d3.1.1] == "T") {
						guideList.split[d3.1.1] = "A"
					}#end if
				}#end d3.1.1 for-loop
				guideList.split
			
				guideList.reverseComplement = "AAAC"

				for (d3.1.2 in length(guideList.split):1) {
					guideList.reverseComplement = paste(guideList.reverseComplement, guideList.split[d3.1.2], sep = "")
				}#end d3.1.2 for-loop
				guideList.reverseComplement
				
				if(d3.1 == 1){
					reverseComplement.list = guideList.reverseComplement
				} else {
					reverseComplement.list = c(reverseComplement.list, guideList.reverseComplement)
				}#end if-else			
				guideList.reverseComplement; guideList; reverseComplement.list
				
			}#end d3.1 for

			guideList = paste("CACC", guideList, sep = ""); guideList
			
			holderMatrix = matrix(ncol = 7, data = c(genebank.matrix[,1], genebank.matrix[,2], genebank.matrix[,5], guideList, reverseComplement.list, genebank.matrix[,4], genebank.matrix[,3]), byrow = FALSE)
			
			if(name.previous != name) {
				compiled.matrix = holderMatrix; #compiled.matrix
				name.previous = name; #name; #name.previous		
			} else {
				compiled.matrix = rbind(compiled.matrix, holderMatrix); #compiled.matrix		
			}#end if-else
		}#end if
	}#end d3 for
	
	for(d5 in 1:length(compiled.matrix[,7])){
		if(grepl(compiled.matrix[d5,7], pattern = "forward") == TRUE){
			compiled.matrix[d5,7] = "+"
		} else if(grepl(compiled.matrix[d5,7], pattern = "reverse") == TRUE){
			compiled.matrix[d5,7] = "-"
		}#end if-else
	}#end d5 for
	
	#saves file to the filteredGuides workFolders[10]
	write.csv(compiled.matrix, file.path(path, workFolders[1], workFolders[10], paste(name, "_filtered_guides.csv", sep = "")), row.names = FALSE, quote = FALSE)	

	files = list.files(file.path(path, workFolders[1], workFolders[10]), full.names = TRUE, pattern = ".csv"); #files

	for (d4 in 1:length(files)){
		guides = read.csv(files[d4]); #guides

		uniques = unique(guides[,5]); #uniques
		
		if(length(guides[,5]) >= 2){
			for (d4.0 in 1:length(uniques)){
				whichOnes = which(guides[,5] == uniques[d4.0]); #whichOnes
				if(length(whichOnes) != 1){
					whichOnes = whichOnes[-1]; whichOnes
					guides = guides[-c(whichOnes),]; #guides				
				} #end if			
			} #end d4.0 for-loop		
		} #end if
		
		filesName = strsplit(files[d4], "/")[[1]]; #filesName
		filesName = filesName[length(filesName)]; #filesName
		
		#save the file in the same folder filteredGuides workFolders[10]
		write.csv(guides, file.path(path, workFolders[1], workFolders[10], paste("unique_", filesName, sep = "")), row.names = FALSE, quote = FALSE)	
	} #end d4 for-loop
	
	#clear vars
	rm(gb.geneList, d1, gb.gene, check, checkExists, gb.geneMatrix, genes.url, genes.names, genes.jobURL, genes.filepaths, name.previous, compiled.matrix, name, genebank, genebank.matrix, files, filesName, guides, uniques, whichOnes, d2, d3, d4, d4.0, d5, lastUnderscore, holderMatrix, entryMITGenebank)
	gc()
	
	entryMITBlat = TRUE
	} else {
		stop("No .gb files were found")
	}#end if
	
}#end entryMITGenebank if
###################

if(entryMITBlat == TRUE){
	#Take the information in the unique_filteredGuides file, and input that into the UCSC BLAT page, to acquire the location of the guide.
	#ucsc BLAT	
	uniqueFilepaths = list.files(file.path(path, workFolders[1], workFolders[10]), full.names = TRUE, pattern = "unique_"); uniqueFilepaths

	for (i0 in 1:length(uniqueFilepaths)){
		
		#using MIT
		uniqueGuidesCSV = read.csv(uniqueFilepaths[i0], header = TRUE); #uniqueGuidesCSV
		
		if (i0 == 1){
			uniqueGuideMatrix = uniqueGuidesCSV
		} else {
			uniqueGuideMatrix = rbind(uniqueGuideMatrix, uniqueGuidesCSV)
		}#end if-else
		
	}#end i0 for
	
	if(file.exists(file.path(path, workFolders[1], workFolders[11], "tempGuidesPostBlat.csv"))){
		blat.resultList = read.csv(file.path(path, workFolders[1], workFolders[11], "tempGuidesPostBlat.csv"), header = TRUE)
		i1 = length(blat.resultList)
	} else {
		i1 = 1
	}#end if-else

	RSDRIVER = rsDriver(port = 4444L, browser = "chrome", check = TRUE, verbose = FALSE)
	RS.browser = RSDRIVER$client

	for (i1 in i1:length(uniqueGuideList)){
	
		RS.browser$navigate("https://genome.ucsc.edu/cgi-bin/hgBlat?command=start")
		Sys.sleep(1)
		
		#css selector for the User Sequence in UCSC BLAT
		#"#firstSection > table:nth-child(1) > tbody:nth-child(1) > tr:nth-child(1) > td:nth-child(1) > table:nth-child(1) > tbody:nth-child(1) > tr:nth-child(1) > td:nth-child(1) > table:nth-child(2) > tbody:nth-child(1) > tr:nth-child(2) > td:nth-child(2) > form:nth-child(4) > table:nth-child(4) > tbody:nth-child(1) > tr:nth-child(3) > td:nth-child(1) > textarea:nth-child(1)"

		ucsc.input = RS.browser$findElement(using = 'css selector', "#firstSection > table:nth-child(1) > tbody:nth-child(1) > tr:nth-child(1) > td:nth-child(1) > table:nth-child(1) > tbody:nth-child(1) > tr:nth-child(1) > td:nth-child(1) > table:nth-child(2) > tbody:nth-child(1) > tr:nth-child(2) > td:nth-child(2) > form:nth-child(4) > table:nth-child(4) > tbody:nth-child(1) > tr:nth-child(3) > td:nth-child(1) > textarea:nth-child(1)")
		ucsc.input$clickElement()
		ucsc.input = RS.browser$sendKeysToActiveElement(list(uniqueGuideMatrix[i1,3]))

		#submit button in UCSC BLAT
		#"#firstSection > table:nth-child(1) > tbody:nth-child(1) > tr:nth-child(1) > td:nth-child(1) > table:nth-child(1) > tbody:nth-child(1) > tr:nth-child(1) > td:nth-child(1) > table:nth-child(2) > tbody:nth-child(1) > tr:nth-child(2) > td:nth-child(2) > form:nth-child(4) > table:nth-child(4) > tbody:nth-child(1) > tr:nth-child(4) > td:nth-child(1) > input:nth-child(1)"

		ucsc.submit = RS.browser$findElement(using = 'css selector', "#firstSection > table:nth-child(1) > tbody:nth-child(1) > tr:nth-child(1) > td:nth-child(1) > table:nth-child(1) > tbody:nth-child(1) > tr:nth-child(1) > td:nth-child(1) > table:nth-child(2) > tbody:nth-child(1) > tr:nth-child(2) > td:nth-child(2) > form:nth-child(4) > table:nth-child(4) > tbody:nth-child(1) > tr:nth-child(4) > td:nth-child(1) > input:nth-child(1)")
		ucsc.submit$clickElement()
		
		Sys.sleep(1)
		
		blat.result = RS.browser$getPageSource()[[1]]
		write(blat.result, file.path(path, workFolders[1], workFolders[11], "BLAT_result.txt"))
		
		Sys.sleep(1)
		
		blat.result = list.files(file.path(path, workFolders[1], workFolders[11]), pattern = "BLAT_result.txt", full.names = TRUE); #blat.result
		blat.result = scan(blat.result, what = "", sep = "\n", quiet = TRUE); #blat.result
		
		check = FALSE		
		while (check == FALSE){
			if (substr(blat.result[1], 1, 5) == " <h2>"){
				check = TRUE
			} else {
				blat.result = blat.result[-1]
			}#end if-else
		}#end while
		
		blat.result = blat.result[15]; #blat.result
		blat.result = strsplit(substr(blat.result, start = (nchar(blat.result)-75), stop = nchar(blat.result)), " ")[[1]]; #blat.result
		blat.result = blat.result[-c(which(blat.result == ""))]; #blat.result
		blat.result = paste("chr", blat.result[6], ":", blat.result[8], "-", blat.result[9], sep = ""); #blat.result
		
		if (i1 == 1){
			blat.resultList = blat.result
		} else {
			blat.resultList = c(blat.resultList, blat.result)
		}#end if-else
		
		#blat.resultList
		
		write.csv(blat.resultList, file.path(path, workFolders[1], workFolders[11], "tempGuidesPostBlat.csv"), row.names = FALSE, quote = FALSE)
		
		Sys.sleep(1)
	}#end i1 for-loop	
	
	uniqueGuideMatrix[,6] = blat.resultList
	
	#ok to remove this file now
	file.remove(file.path(path, workFolders[1], workFolders[11], "tempGuidesPostBlat.csv")); file.remove(file.path(path, workFolders[1], workFolders[11], "BLAT_result.txt"))
	
	#save file to BLAT workFolders[11]
	write.csv(uniqueGuideMatrix, file = file.path(path, workFolders[1], workFolders[11], "guidesPostBLAT.csv"), row.names = FALSE, quote = FALSE)
	
	#close the browser
	RS.browser$close()
	RSDRIVER$server$stop()
	
	#clear vars
	rm(uniqueFilepaths, i0, i1, uniqueGuidesCSV, ucsc.input, ucsc.submit, blat.result, check, blat.resultList, uniqueGuideMatrix, entryMITBlat); gc()
	
	entryPostBlat = TRUE
}#end entryMITBlat if

}#END MIT MODULE 
##########

##########
#CCTOP
if(cctop == TRUE & mit.edu == FALSE){
if(entryCCTOPSubmit == TRUE){
#cctop css selectors

	#batch mode element
		#"#sendQuery > table:nth-child(4) > tbody:nth-child(1) > tr:nth-child(1) > th:nth-child(2) > input:nth-child(1)"
		
	#name element
		#"#name"
		
	#fasta file element
		#"#fastaFile"

	#email element
		#"#email"
		
	#U6 element
		#"#sendQuery > table:nth-child(8) > tbody:nth-child(1) > tr:nth-child(2) > td:nth-child(1) > div:nth-child(10) > input:nth-child(2)"
	
	#species pull-down element
		#"#species > option:nth-child(19)"
	#zip file download element
		#"#middleColumn > a:nth-child(3)"	
	
	#list the compiledFa.fa files in the compiledFastas workFolders[4]
	files.compiledFa = list.files(file.path(path, workFolders[1], workFolders[4]), full.names = TRUE, pattern = "compiled.fa")

	if(length(files.compiledFa) >= 1){
	#open the browser
	RSDRIVER = rsDriver(port = 4444L, browser = "chrome", check = TRUE, verbose = FALSE)
	RS.browser = RSDRIVER$client
	
	for (f0 in 1:length(files.compiledFa)){	
		#navigate to CCTOP
		RS.browser$navigate("http://crispr.cos.uni-heidelberg.de/")

		#select batch mode
		cctop.batchMode = RS.browser$findElement(using = 'css selector', "#sendQuery > table:nth-child(4) > tbody:nth-child(1) > tr:nth-child(1) > th:nth-child(2) > input:nth-child(1)")
		cctop.batchMode$clickElement()
		
		#enter a name into the field
		cctop.nameToSend = strsplit(files.compiledFa[f0], "/")[[1]]; cctop.nameToSend
		cctop.nameToSend = cctop.nameToSend[length(cctop.nameToSend)]; cctop.nameToSend
		cctop.nameToSend = strsplit(cctop.nameToSend, "_")[[1]]
		cctop.nameToSend = paste(cctop.nameToSend[1], cctop.nameToSend[2], sep = "_")
		#cctop.nameToSend = substr(cctop.nameToSend, start = 1, stop = (nchar(cctop.nameToSend)-14)); cctop.nameToSend
		cctop.name = RS.browser$findElement(using = 'css selector', "#name")
		cctop.name$clearElement()
		cctop.name$sendKeysToElement(list(cctop.nameToSend))		
		
		#send the file to the openFileDialog element
		cctop.fastaFile = RS.browser$findElement(using = 'css selector', "#fastaFile")
			tempFasta = tempfile(fileext = ".fa")
			fastaFile = scan(file = files.compiledFa[f0], what = "", sep = "\n", quiet = TRUE)
			write(fastaFile, tempFasta)
		cctop.fastaFile$sendKeysToElement(list(tempFasta))

		#enter the email address
		cctop.email = RS.browser$findElement(using = 'css selector', "#email")
		cctop.email$clearElement()
		cctop.email$sendKeysToActiveElement(list(email))
		
		#set promoter to U6
		cctop.promotor = RS.browser$findElement(using = 'css selector', "#sendQuery > table:nth-child(8) > tbody:nth-child(1) > tr:nth-child(2) > td:nth-child(1) > div:nth-child(10) > input:nth-child(2)")
		cctop.promotor$clickElement()
		
		#set pull-down menu to Human genome (Hg38)
		cctop.pullDown = RS.browser$findElement(using = 'css selector', "#species > option:nth-child(19)")
		cctop.pullDown$clickElement()
		

		#click submit
		cctop.submit = RS.browser$findElement(using = 'css selector', "#sendQuery > input:nth-child(11)")
		cctop.submit$clickElement()
	
		#get the page URL
		Sys.sleep(3)
		if (f0 == 1){
			cctop.url = RS.browser$getCurrentUrl()[[1]]; #cctop.url
		} else {
			cctop.url = c(cctop.url, RS.browser$getCurrentUrl()[[1]])
		} #end if-else	
		
	} #end f0 for		
	
	#close the browser and stop the server
	RS.browser$close()
	RSDRIVER$server$stop()
	
	#save to the CCTOP workFolders[7]
	write(cctop.url, file.path(path, workFolders[1], workFolders[7], "cctop_jobSatusURLS.txt"))
	
	#clear vars
	rm(files.compiledFa, f0, cctop.batchMode, cctop.nameToSend, cctop.name, cctop.fastaFile, tempFasta, fastaFile, cctop.email, cctop.promotor, cctop.pullDown, cctop.submit, cctop.url, entryCCTOPSubmit)
	gc()
	
	entryCCTOPMonitor = TRUE
	} else {
		stop("No compiled.fa files were found")
	}#end if
}#end entryCCTOPSubmit if
#############

###############		
		
#check the job status page by refreshing every so often, and completion can be measured by comparing the two urls, past and present. Once complete, download the zip file.
if(entryCCTOPMonitor == TRUE){
	cctop.url = list.files(file.path(path, workFolders[1], workFolders[7]), full.names = TRUE, pattern = "cctop_jobSatusURLS.txt")
	
	if(length(cctop.url) >= 1){
	cctop.urlFile = scan(cctop.url, what = "", sep = "\n", quiet = TRUE)
	
	#open the browser
	RSDRIVER = rsDriver(port = 4444L, browser = "chrome", check = TRUE, verbose = FALSE)
	RS.browser = RSDRIVER$client
	
	for(g0 in 1:length(cctop.urlFile)){
		RS.browser$navigate(cctop.urlFile[g0])

		check = FALSE
		while (check == FALSE){
			#get the text from the status text on the page
				#for checking for the download
			#css selector for text
			# "#middleColumn > h1:nth-child(1)"
			
			cctop.statusText = RS.browser$findElement(using = 'css', "#middleColumn > h1:nth-child(1)")
			cctop.statusText = cctop.statusText$getElementText()[[1]]; cctop.statusText
			if (length(grep(cctop.statusText, pattern = "Results")) == 1){
				check = TRUE
				cctop.zipName = strsplit(cctop.statusText, " ")[[1]]; cctop.zipName
				cctop.zipName = cctop.zipName[length(cctop.zipName)]
				
				#detect url for zip file	
				cctop.zip = RS.browser$findElement(using = 'css', "#middleColumn > a:nth-child(3)")
				cctop.zip = cctop.zip$getElementAttribute('href')[[1]]; #cctop.zip
				
				Sys.sleep(1)
				
				#download the zip file to the zip workFolders[8], unless it already exists
				if(!file.exists(file.path(path, workFolders[1], workFolders[7], workFolders[8], paste("CCTOP_Results_", cctop.zipName, ".zip", sep = "")))){
					download.file(url = cctop.zip, destfile = file.path(path, workFolders[1], workFolders[7], workFolders[8], paste("CCTOP_Results_", cctop.zipName, ".zip", sep = "")))
					
					print(paste("Downloaded ", cctop.zipName, ".", sep = ""))
				} else {
					print(paste("CCTOP_Results_", cctop.zipName, ".zip", " already exists.", sep = ""))				
				}#end if-else
			} else {
				print(paste(Sys.time(), "  :  Nothing just yet...", sep = ""))
				Sys.sleep(178)
				RS.browser$refresh()
				Sys.sleep(1)
			}#end if
		}#end while
	} #end g0 for
	
	#close the browser and stop the server
	RS.browser$close()
	RSDRIVER$server$stop()
	
	rm(cctop.url, cctop.urlFile, g0, cctop.statusText, cctop.zip, cctop.zipName)
	gc()
##########

#CCTOP
	#extract files from a zip to the XLS workFolders[9]
	extractDir = file.path(path, workFolders[1], workFolders[7], workFolders[9]); #extractDir
	extractFile = list.files(file.path(path, workFolders[1], workFolders[7], workFolders[8]), full.names = TRUE, pattern = ".zip"); #extractFile
	
	for (g1 in 1: length(extractFile)){
		unzip(zipfile = extractFile[g1], exdir = extractDir)
	} #end g1 for

	rm(extractDir, extractFile, g1)

###############	

	rm(entryCCTOPMonitor)
	
	entryCCTOP.XLS = TRUE
	
	} else {
		stop("No compiled.fa files were found")
	}#end if
	
}#end entryCCTOPMonitor if	

##########
#cctop score calculator
#uses the Zhang score calclation algorithm
if (entryCCTOP.XLS == TRUE){
	#weighted value set
	W.e = c(0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583)
	invW.e = 1 - W.e
	
	#picks up the csv files in the CCTOP XLS workFolders[9]
	XLS.list = list.files(file.path(path, workFolders[1], workFolders[7], workFolders[9]), full.names = TRUE, pattern = ".xls", recursive = TRUE); #XLS.list
	#XLS.list = list.files(file.path(path, workFolders[1], workFolders[7]), full.names = TRUE, pattern = ".csv"); #XLS.list
	
	if(length(XLS.list) >= 1){
	for (h0 in 1:length(XLS.list)){
		
		#acquire gene name
		XLS.name = strsplit(XLS.list[h0], "/")[[1]]
		XLS.name = XLS.name[length(XLS.name)]
		XLS.name = substr(XLS.name, start = 1, stop = (nchar(XLS.name)-4)); #XLS.name	
		
		#groom the file
		XLS = scan(as.character(XLS.list[h0]), what = "", sep = "\n", quiet = TRUE)
		XLS.input = substr(XLS[2], 8, nchar(XLS[2]))
		XLS = XLS[-c(1:9)]; #head(XLS, 50)
		XLS = XLS[-c(grep(XLS, pattern = "Oligo substituting"))]; #length(XLS)
		
		XLS.targets = which(substr(XLS, 1, 1) == "T"); #XLS.targets
		
		for (h0.0 in XLS.targets) {
			targetSequence = strsplit(XLS[h0.0], split = "\t")[[1]][2]; #targetSequence
			oligoFwd = strsplit(XLS[h0.0 + 1], split = "\t")[[1]][2]; #oligoFwd
			oligoRev = strsplit(XLS[h0.0 + 2], split = "\t")[[1]][2]; #oligoRev
			
			#accounts for the last target in the list
			if (h0.0 != XLS.targets[length(XLS.targets)]){
				offTargetList = XLS[c((h0.0 + 4):(XLS.targets[(which(XLS.targets == h0.0) + 1)] - 1))]; #print(offTargetList)
			} else {
				offTargetList = XLS[c((h0.0 + 4):length(XLS))]; #print(offTargetList)
			} #end if-else
			
			#actual score calculator
			for (h0.1.1 in 1:length(offTargetList)){
				if (h0.1.1 == 1) {
					seqAndAlignment = strsplit(offTargetList[h0.1.1], split = "\t")[[1]][c(6,8)]; #print(seqAndAlignment)
					guidePosition =  strsplit(offTargetList[h0.1.1], split = "\t")[[1]][c(1:4)]; #print(guidePosition)
					
					if(guidePosition[4] == "+"){
						guidePosition[3] = as.integer(guidePosition[3]) - 3			
					} else if (guidePosition[4] == "-") {
						guidePosition[2] = as.integer(guidePosition[2]) + 3			
					} #end if-else
					
					genomePosition = paste(guidePosition[1], ":", guidePosition[2], "-", guidePosition[3], sep = ""); #print(genomePosition)
					strandedness = guidePosition[4]; #print(strandedness)
					
					target = seqAndAlignment[2]; #print(target)
					
					target = strsplit(target, "")[[1]]; target
					target = target[-c(which(target == "["), 22:25)]; target

					summedSingleGuideScores = 100
				} else {
					seqAndAlignment = strsplit(offTargetList[h0.1.1], split = "\t")[[1]][c(6,8)]; #print(seqAndAlignment)
					
					offTarget = seqAndAlignment[2]; #print(offTarget)

					offTarget = strsplit(offTarget, "")[[1]]; offTarget
					offTarget = offTarget[-c(which(offTarget == "["), 22:25)]; offTarget

					e = which(offTarget[] != target[]); #print(e)
					nmm = length(e); #print(nmm)

					if (nmm > 1){ #number of mismatches is more than 1
						d = ((max(e) - min(e))/(nmm - 1))
					} else { #number of mismatches is less than or equal to 1
						d = 0
					} #end if-else
					#print(d)

					productW.e = 1

					for (h0.1.1.1 in 1:nmm){
						productW.e = productW.e * invW.e[e[h0.1.1.1]]
					} #end h0.1.1.1 for
					#print(productW.e); productW.e
					
					#catches the event that the off target is identical to the guide
					if (identical(productW.e, numeric(0))){
						singleGuideScore = 100
					} else {						
						singleGuideScore = 100 * productW.e * (1/((4*((19 - d)/19)) + 1)) * (1/(nmm^2))
					}#end if-else					
					
					summedSingleGuideScores = summedSingleGuideScores + singleGuideScore
				} #end if-else		
			} #end h0.1.1 for
			
			aggregateScore = 100 * (100 / summedSingleGuideScores); #aggregateScore
			
			if (h0.0 == 1) {
				targetAndScoreData = c(XLS.name, aggregateScore, targetSequence, oligoFwd, oligoRev, genomePosition, strandedness)
			} else {
				targetAndScoreData = c(targetAndScoreData, XLS.name, aggregateScore, targetSequence, oligoFwd, oligoRev, genomePosition, strandedness)
			} #end if-else
			
		} #end h0.1 for
		
		targetAndScoreMatrix = matrix(ncol = 7, data = targetAndScoreData, byrow = TRUE); targetAndScoreMatrix
		
		#save targetAndScoreMatrix to the XLS workFolders[9]
		write.csv(targetAndScoreMatrix, file.path(path, workFolders[1], workFolders[7], workFolders[9], paste("CCTOP_Results_", XLS.name, ".csv", sep = "")), row.names = FALSE, quote = FALSE)
		
	}#end h0 for	
	
	rm(W.e, invW.e, XLS, XLS.input, XLS.targets, targetSequence, oligoFwd, oligoRev, offTargetList, seqAndAlignment, guidePosition, genomePosition, strandedness, target, summedSingleGuideScores, offTarget, e, nmm, d, productW.e, singleGuideScore, aggregateScore, targetAndScoreMatrix, targetAndScoreData, h0, h0.1.1, h0.1.1.1)
	
	} else {
		stop("No .xls files were found")
	}#end if

######
#filter the guides from the CCTOP results that were just output.
	cctopResults.filePath = list.files(file.path(path, workFolders[1], workFolders[7], workFolders[9]), pattern = "CCTOP_Results_", full.names = TRUE)

	#dummy vectors
	name.previous = ""; #name.previous
	name.next = ""; #name.next
	compiledMatrix = matrix(nrow = 1, ncol = 7); #compiledMatrix

	for (h1 in 1:length(cctopResults.filePath)){		
		targetAndScoreMatrix = read.csv(cctopResults.filePath[h1], header = TRUE, colClasses = c('character', 'numeric', rep('character', 5))); #targetAndScoreMatrix
		
		#acquire gene name
		name = strsplit(cctopResults.filePath[h1], split = "/")[[1]]; #name
		name = name[length(name)]; #name
		name = strsplit(name, split = "_")[[1]][3]; #name
		
		if(h1 != length(cctopResults.filePath)){
			name.next = strsplit(cctopResults.filePath[h1+1], split = "/")[[1]]; #name.next
			name.next = name.next[length(name.next)]; #name.next
			name.next = strsplit(name.next, split = "_")[[1]][3]; #name.next
		} else {
			#if the next name in the list doesn't exist bc it's the end of the list, make name.next empty
			name.next = ""
		}#end if
		
		#filter guides according to threshold value
		targetAndScoreMatrix = targetAndScoreMatrix[-c(which(targetAndScoreMatrix[,2] < threshold)),]; #targetAndScoreMatrix
		
		#keep same gene guides together
		if(name.previous != name) {
			compiledMatrix = targetAndScoreMatrix
			name.previous = name; #name; #name.previous		
		} else {
			compiledMatrix = rbind(compiledMatrix, targetAndScoreMatrix); #compiled.matrix		
		}#end if-else

		if(name.next != name){
			#saves file to the filteredGuides workFolders[10]
			write.csv(compiledMatrix, file.path(path, workFolders[1], workFolders[10], paste(name, "_filtered_guides.csv", sep = "")), row.names = FALSE, quote = FALSE)
		}#end if		
	}#end h1 for
	
	#clear vars
	rm(cctopResults.filePath, name.previous, name.next, compiledMatrix, guidesPostBLAT, h1, targetAndScoreMatrix, name)
#################

######aggregate the unique filtered guides into guidesPostBLAT.csv and unique_filteredGuides.csv
	filteredGuides.filePath = list.files(file.path(path, workFolders[1], workFolders[10]), pattern = "filtered_guides.csv", full.names = TRUE); filteredGuides.filePath
	
	for(h2 in 1:length(filteredGuides.filePath)){
		filteredGuidesMatrix = read.csv(filteredGuides.filePath[h2], header = TRUE, colClasses = c('character', 'numeric', rep('character', 5))); #head(filteredGuidesMatrix, 20)	

		uniques = unique(filteredGuidesMatrix[,5]); #uniques
		
		if(length(filteredGuidesMatrix[,5]) >= 2){
			for (h2.0 in 1:length(uniques)){
				whichOnes = which(filteredGuidesMatrix[,5] == uniques[h2.0]); #whichOnes
				if(length(whichOnes) != 1){
					whichOnes = whichOnes[-1]; #whichOnes
					filteredGuidesMatrix = filteredGuidesMatrix[-c(whichOnes),]; #length(filteredGuidesMatrix[,5])				
				} #end if			
			} #end h2.0 for		
		} #end if
		
		filesName = strsplit(filteredGuides.filePath[h2], "/")[[1]]; #filesName
		filesName = filesName[length(filesName)]; #filesName
		
		#get rid of repeat guides due to overlapping of files
		for(h2.1 in 1:length(filteredGuidesMatrix[,3])){
			guidesRepeats = which(filteredGuidesMatrix[,3] == filteredGuidesMatrix[h2.1,3]); #guidesRepeats
			
			if(length(guidesRepeats) > 1){
				guidesRepeats = guidesRepeats[-1]
				filteredGuidesMatrix = filteredGuidesMatrix[-c(guidesRepeats),]
			}#end if
		}#end h2.1 for
		#print(length(filteredGuidesMatrix[,3]))
		
		#save the file in the same folder filteredGuides workFolders[10]
		write.csv(filteredGuidesMatrix, file.path(path, workFolders[1], workFolders[10], paste("unique_", filesName, sep = "")), row.names = FALSE, quote = FALSE)		
	}#end h2 for	
	rm(guidesRepeats)
	
	unique_filteredGuides.filePath = list.files(file.path(path, workFolders[1], workFolders[10]), pattern = "unique_", full.names = TRUE)
	
	for(h3 in 1: length(unique_filteredGuides.filePath)){
		unique_filteredGuidesMatrix = read.csv(unique_filteredGuides.filePath[h3], header = TRUE, colClasses = c('character', 'numeric', rep('character', 5))); #head(unique_filteredGuidesMatrix, 50)

		if(h3 == 1){
			guidesPostBLAT = unique_filteredGuidesMatrix; #head(guidesPostBLAT, 20)
		} else {
			guidesPostBLAT = rbind(guidesPostBLAT, unique_filteredGuidesMatrix)		
		}#end if-else
	}#end h3 for
	
	#get rid of repeat guides due to overlapping of files, or any other non-unique ones that snuck in
	for(h3.0 in 1:length(guidesPostBLAT[,3])){
		guidesRepeats = which(guidesPostBLAT[,3] == guidesPostBLAT[h3.0,3]); #print(guidesRepeats)
		
		if(length(guidesRepeats) > 1){
			guidesRepeats = guidesRepeats[-1]
			guidesPostBLAT = guidesPostBLAT[-c(guidesRepeats),]
		}#end if
	}#end h3.0 for
	#print(length(guidesPostBLAT[,3]))

	
	#save to same filteredGuides workFolders[10]
	write.csv(guidesPostBLAT, file.path(path, workFolders[1], workFolders[11], "guidesPostBLAT.csv"), row.names = FALSE, quote = FALSE)		
	
	#clear vars
	rm(filteredGuides.filePath, h2, filteredGuidesMatrix, uniques, h2.0, whichOnes, filesName, unique_filteredGuides.filePath, guidesPostBLAT, h2.1, unique_filteredGuidesMatrix, h3, guidesRepeats, entryCCTOP.XLS)
	gc()
	
	entryPostBlat = TRUE
	
}#end entryCCTOP.XLS if

}#END CCTOP MODULE		
#############	

##############
#take the location, input that into the UCSC Get DNA page, and acquire the sequence around the guide.
#copy the ucsc get dna code, and change the padding to something like 250 on either side.
if(entryPostBlat == TRUE){	
	j0 = 1
	guidesPostBLAT = list.files(file.path(path, workFolders[1], workFolders[11]), pattern = "guidesPostBLAT.csv", full.names = TRUE); #guidesPostBLAT
	
	padded250bpList = list.files(file.path(path, workFolders[1], workFolders[11]), pattern = "padded250bp.fa", full.names = TRUE); #head(padded250bpList, 20)

	#figures out how many guide fastas have already been downloaded
	if(length(padded250bpList) > 0){ 
		for(j2 in 1:length(padded250bpList)){
			padded250bpNumber = strsplit(padded250bpList[j2], split = "/")[[1]]; #padded250bpNumber
			padded250bpNumber = padded250bpNumber[length(padded250bpNumber)]; #padded250bpNumber
			padded250bpNumber = strsplit(substr(padded250bpNumber, start = 1, stop = 7), split = "")[[1]]; #padded250bpNumber
			
			checkNumber = FALSE
			while(checkNumber == FALSE){
				if(padded250bpNumber[1] == "0"){
					padded250bpNumber = padded250bpNumber[-1]; #print(padded250bpNumber)
				} else {
					checkNumber = TRUE
				}#end if-else
			}#end while
			
			padded250bpNumber.0 = ""
			for(j2.0 in 1:length(padded250bpNumber)){
				#handles numbers less than 10
				if(length(padded250bpNumber) != 1){
						padded250bpNumber.0 = paste(padded250bpNumber.0, padded250bpNumber[j2.0], sep = "")
				} else {
				padded250bpNumber.0 = padded250bpNumber
				}#end if-else
			}#end j2.0 for
			
			#padded250bpNumber.0
			
			j0 = as.integer(padded250bpNumber.0) + 1; #j0
		}#end j2 for
		
		print(paste("Starting at #", j0, " in guidesPostBLAT.csv", sep = ""))
	}#end if
	
	if(length(guidesPostBLAT) == 1){
	guidesPostBLAT = read.csv(guidesPostBLAT, header = TRUE, colClasses = c(rep('character', 3))); #guidesPostBLAT
	
	#print estimated time to completion
	print(paste("Estimated ", round((length(guidesPostBLAT[,1])-j0)*0.00291666666, 2), " hours to complete ", length(guidesPostBLAT[,1]) - j0, " of ", length(guidesPostBLAT[,1]), " FASTAs for primer3", sep = ""))
	
	#link to the ucsc gateway
	# "http://genome.ucsc.edu/cgi-bin/hgGateway"
	
	RSDRIVER = rsDriver(port = 4444L, browser = "chrome", check = TRUE, verbose = FALSE)
	RS.browser = RSDRIVER$client

	for (j0 in j0:length(guidesPostBLAT[,2])){
		RS.browser$navigate("http://genome.ucsc.edu/cgi-bin/hgGateway")
		
		Sys.sleep(1)

		#css selector for position input for ucsc
		#"#positionInput"
		ucsc.input = RS.browser$findElement(using = 'css', "#positionInput")
		ucsc.input$clickElement()
		ucsc.input = RS.browser$sendKeysToActiveElement(list(guidesPostBLAT[j0,3]))
		
		#css selector for submit button
		#".jwGoButton"
		ucsc.submit = RS.browser$findElement(using = 'css', '.jwGoButton')
		ucsc.submit$clickElement()
		
		#css selector for viewDNA
		#"#dnaLink"
		ucsc.viewDNA = RS.browser$findElement(using = 'css', "#dnaLink")
		RS.browser$navigate(ucsc.viewDNA$getElementAttribute(attr = 'href')[[1]])
		
		#xpath selector for 5' padding in getDNA
		# '//*[@id="hgSeq.padding5"]'
		ucsc.getDNA5 = RS.browser$findElement(using = 'xpath', '//*[@id="hgSeq.padding5"]')
		ucsc.getDNA5$setElementAttribute(attributeName = 'value', 250)
		
		#xpath for 3' padding in getDNA
		#'//*[@id="hgSeq.padding3"]'
		ucsc.getDNA3 = RS.browser$findElement(using = 'xpath', '//*[@id="hgSeq.padding3"]')
		ucsc.getDNA3$setElementAttribute(attributeName = 'value', 250)
		
		#css for submit
		#"#firstSection > table:nth-child(1) > tbody:nth-child(1) > tr:nth-child(1) > td:nth-child(1) > table:nth-child(1) > tbody:nth-child(1) > tr:nth-child(1) > td:nth-child(1) > table:nth-child(2) > tbody:nth-child(1) > tr:nth-child(2) > td:nth-child(2) > form:nth-child(5) > p:nth-child(30) > input:nth-child(1)"
		ucsc.submit = RS.browser$findElement(using = 'css', "#firstSection > table:nth-child(1) > tbody:nth-child(1) > tr:nth-child(1) > td:nth-child(1) > table:nth-child(1) > tbody:nth-child(1) > tr:nth-child(1) > td:nth-child(1) > table:nth-child(2) > tbody:nth-child(1) > tr:nth-child(2) > td:nth-child(2) > form:nth-child(5) > p:nth-child(30) > input:nth-child(1)")
		ucsc.submit$clickElement()
		
		#zero padding module for keeping guides in order
		if (j0 < 10){ #ones
			zeroPad = "000000"
			guideNumber = paste(zeroPad, j0, sep = "")
		} else if (j0 >= 10 & j0 < 100) { #tens
			zeroPad = "00000"
			guideNumber = paste(zeroPad, j0, sep = "")
		} else if (j0 >= 100 & j0 < 1000) { #thousands
			zeroPad = "0000"
			guideNumber = paste(zeroPad, j0, sep = "")
		} else if (j0 >= 1000 & j0 < 10000) { #ten thousands
			zeroPad = "000"
			guideNumber = paste(zeroPad, j0, sep = "")
		} else if (j0 >= 10000 & j0 < 100000) { #hundred thousands
			zeroPad = "00"
			guideNumber = paste(zeroPad, j0, sep = "") 
		} else if (j0 >= 100000 & j0 < 1000000) { #millions
			zeroPad = "0"
			guideNumber = paste(zeroPad, j0, sep = "")
		} else if (j0 >= 1000000 & j0 < 10000000) { #ten millions
			zeroPad = ""
			guideNumber = paste(zeroPad, j0, sep = "")		
		}#end zero padding if-else
		
		#print(guideNumber)
		
		download.file(RS.browser$getCurrentUrl()[[1]], file.path(path, workFolders[1], workFolders[11],	paste(guideNumber, "_", guidesPostBLAT[j0,1], "_padded250bp.fa", sep = "")))
		
		Sys.sleep(1)
	} #end j0 for
	
	RS.browser$close()
	RSDRIVER$server$stop()
	
	ucsc.viewDNA.filePath =  list.files(file.path(path, workFolders[1], workFolders[11]), full.names = TRUE, pattern = "padded250bp.fa"); #ucsc.viewDNA.filePath
	
	for (j1 in 1:length(ucsc.viewDNA.filePath)){
		ucsc.viewDNA.file = scan(ucsc.viewDNA.filePath[j1], what = "", sep = "\n", quiet = TRUE)
		
		check = FALSE
		
		while (check == FALSE) {
			if(substr(ucsc.viewDNA.file[1], 1, 5) == ">hg38") {
				check = TRUE
			}#end if	
			else {
				ucsc.viewDNA.file = ucsc.viewDNA.file[-1]
			}#end else
		}#end while-loop
		
		ucsc.viewDNA.fileName = strsplit(ucsc.viewDNA.filePath[j1], "/")[[1]]; ucsc.viewDNA.fileName
		ucsc.viewDNA.fileName = ucsc.viewDNA.fileName[length(ucsc.viewDNA.fileName)]; ucsc.viewDNA.fileName
		ucsc.viewDNA.fileName = substr(ucsc.viewDNA.fileName, start = 1, stop = (nchar(ucsc.viewDNA.fileName)-3)); ucsc.viewDNA.fileName
		ucsc.viewDNA.file[1] = paste(">", ucsc.viewDNA.fileName, sep = ""); ucsc.viewDNA.file
		
		ucsc.viewDNA.file = ucsc.viewDNA.file[-c(which(ucsc.viewDNA.file == "</PRE>"):length(ucsc.viewDNA.file))]
		write(ucsc.viewDNA.file, ucsc.viewDNA.filePath[j1])
	} #end j1 for-loop

	#clear vars
	rm(guidesPostBLAT, j0, ucsc.input, ucsc.submit, ucsc.viewDNA, ucsc.getDNA5, ucsc.getDNA3, zeroPad, guideNumber, ucsc.viewDNA.filePath, j1, ucsc.viewDNA.file, check, ucsc.viewDNA.fileName, entryPostBlat, padded250bpList, j2, padded250bpNumber, padded250bpNumber.0, j2.0)
	gc()
	
	entryPrimer3 = TRUE
	
	} else {
		stop("No guidesPostBLAT.csv file was found")
	}#end if

}#end entryPostBlat if
################

################
#take the sequence that was acquired from Get DNA, and input it into the Primer3 website (http://primer3.ut.ee/)
#get the primers from Primer3, and deposit it in a file with the highest scoring guides.
if(entryPrimer3 == TRUE){
	primer3.filePath = list.files(file.path(path, workFolders[1], workFolders[11]), pattern = "padded250bp", full.names = TRUE); length(primer3.filePath)
	
	if(length(primer3.filePath) >= 1){
		k0 = 1
		
		if(file.exists(file.path(path, workFolders[1], workFolders[12], "primerMatrix.csv"))){
			primer3.matrix = read.csv(file.path(path, workFolders[1], workFolders[12], "primerMatrix.csv"), header = TRUE, colClasses = c(rep('character', 4))); #print(primer3.matrix)
			k0 = length(primer3.matrix[,1])+1	
		}#end if	
			
		print(paste("Starting at #", k0, " in primerMatrix.csv", sep = ""))
		
		RSDRIVER = rsDriver(port = 4444L, browser = "chrome", check = TRUE, verbose = FALSE)
		RS.browser = RSDRIVER$client
		
		print(paste("Estimated ", round((length(primer3.filePath)-k0)*0.008333333, 2), "hours to complete primerMatrix.csv", sep = ""))
		
		for (k0 in k0:length(primer3.filePath)){
			RS.browser$navigate("http://bioinfo.ut.ee/primer3-0.4.0/")
			
			primer3.file = scan(primer3.filePath[k0], what = "", sep = "\n", quiet = TRUE); #primer3.file
			
			for (k0.0 in 2:length(primer3.file)){
				if (k0.0 == 2) {
					primer3.file1 = primer3.file[k0.0]
				} else {
					primer3.file1 = paste(primer3.file1, primer3.file[k0.0], sep = "")
				} #end if
			} #end k0.0 for
			#primer3.file
			
			#css selector for textarea
			#"body > form:nth-child(1) > p:nth-child(6) > textarea:nth-child(1)"
			primer3.textarea = RS.browser$findElement(using = 'css', "body > form:nth-child(1) > p:nth-child(6) > textarea:nth-child(1)")
			primer3.textarea$clickElement()
			primer3.textarea = RS.browser$sendKeysToActiveElement(list(primer3.file1))
		
			#css selector for target
			#"body > form:nth-child(1) > table:nth-child(9) > tbody:nth-child(1) > tr:nth-child(2) > td:nth-child(2) > input:nth-child(1)"
			primer3.target = RS.browser$findElement(using = 'css', "body > form:nth-child(1) > table:nth-child(9) > tbody:nth-child(1) > tr:nth-child(2) > td:nth-child(2) > input:nth-child(1)")
			primer3.target$setElementAttribute(attributeName = 'value', '251,23')
			
			#css selector for submit
			#"body > form:nth-child(1) > p:nth-child(12) > input:nth-child(1)"
			primer3.submit = RS.browser$findElement(using = 'css', "body > form:nth-child(1) > p:nth-child(12) > input:nth-child(1)")
			primer3.submit$clickElement()
			
			Sys.sleep(1)
			
			#css selector for primer3 output
			#body > pre:nth-child(3) > a:nth-child(8)
			primer3.output = RS.browser$getPageSource()[[1]]
			
			#save file in primer3 workFolders[12]
			write(primer3.output, file.path(path, workFolders[1], workFolders[12], "primer3_output.txt"))		
			primer3.outputFilePath = list.files(file.path(path, workFolders[1], workFolders[12]), pattern = "primer3_output.txt", full.names = TRUE); #primer3.outputFilePath
			
			guideNumber = strsplit(primer3.filePath[k0], "/")[[1]]; #guideNumber
			guideNumber = guideNumber[length(guideNumber)]; #guideNumber
			guideNumber = substr(guideNumber, start = 1, stop = 8); #guideNumber			
			
			primer3.file = scan(primer3.outputFilePath, what = "", sep = "\n", quiet = TRUE); #primer3.file

			if(substr(primer3.file[12], start = 1, stop = 4) == "LEFT"){
				primer3.file = primer3.file[c(12,13,16)]; #primer3.file
				leftPrimer = strsplit(primer3.file[1], " ")[[1]]; #leftPrimer
				leftPrimer = leftPrimer[length(leftPrimer)]; #leftPrimer
				rightPrimer = strsplit(primer3.file[2], " ")[[1]]; #rightPrimer
				rightPrimer = rightPrimer[length(rightPrimer)]; #rightPrimer
				productSize = strsplit(primer3.file[3], " ")[[1]]; #productSize
				productSize = productSize[3]
			} else {
				leftPrimer = NA
				rightPrimer = NA
				productSize = NA
				print(paste(substr(primer3.file[35], start = 5, stop = 25), " for guide #", guideNumber, sep = ""))
			}#end if-else
			
			primer3.result = c(guideNumber, leftPrimer, rightPrimer, productSize)
			
			if (k0 == 1) {
				primer3.matrix = matrix(ncol = 4, data = primer3.result, byrow = TRUE)
			} else {
				primer3.matrix = rbind(primer3.matrix, primer3.result)
			} #end if-else
			
			#save to primer3 workFolders[12]
			write.csv(primer3.matrix, file.path(path, workFolders[1], workFolders[12], "primerMatrix.csv"), row.names = FALSE, quote = FALSE)
			
		} #end k0 for-loop
		
		RS.browser$close()
		RSDRIVER$server$stop()	
		
		#clear vars
		rm(primer3.filePath, k0, primer3.file, k0.0, primer3.file1, primer3.textarea, primer3.target, primer3.submit, primer3.output, leftPrimer, rightPrimer, productSize, primer3.result, primer3.matrix, guideNumber, primer3.outputFilePath, entryPrimer3)
		
		entryResults = TRUE
		
	} else {
		stop("No padded250bp.fa files were found")
	}#end if

}#end entryPrimer3 if	
####################	

#################
#now aggregate everything!
if(entryResults == TRUE){
	primerMatrix.filePath = list.files(file.path(path, workFolders[1], workFolders[12]), pattern = "primerMatrix", full.names = TRUE); primerMatrix.filePath
	
	guidesPostBLAT.filePath = list.files(file.path(path, workFolders[1], workFolders[11]), pattern = "guidesPostBLAT", full.names = TRUE); guidesPostBLAT.filePath
	
	if(length(primerMatrix.filePath) == 1 & length(guidesPostBLAT.filePath) == 1){
		primerMatrix = read.csv(primerMatrix.filePath, header = TRUE); #head(primerMatrix, 50); length(primerMatrix)
			
		guidesPostBLAT = read.csv(guidesPostBLAT.filePath, header = TRUE); #head(guidesPostBLAT, 20)
		
		#report matrix structure, 11 columns
			#gene, guide score, guide sequence + CACCg, reverse complement + CAAA, guide position, guide strandedness, left primer sequence, right primer sequence, PCR product size
		reportMatrix = matrix(ncol = 11, nrow = length(guidesPostBLAT[,1]), data = c(guidesPostBLAT[,1], guidesPostBLAT[,2], guidesPostBLAT[,3], guidesPostBLAT[,4], guidesPostBLAT[,5], guidesPostBLAT[,6], guidesPostBLAT[,7], primerMatrix[,2], primerMatrix[,3], primerMatrix[,4], primerMatrix[,1]), byrow = FALSE); #reportMatrix
		
		write.csv(reportMatrix, file.path(path, workFolders[13], "FINAL RESULTS.csv"), row.names = FALSE, quote = FALSE)

		#clear vars
		rm(primerMatrix.filePath, unique_filteredGuides.filePath, guidesPostBLAT.filePath, primerMatrix, l0, unique_filteredGuides.CSV, unique_filteredGuides, guidesPostBLAT, reportMatrix, entryResults)
		
	} else {
		stop("Necessary input file missing to generate results")
	}#end if
	
}#end entryResults if	

	gc()
	print("Thanks for using gRNA-wrapR! You the real MVP!"); print(paste("Job completed at:  ", Sys.time(), sep = ""))
}#end gRNA.wrapR function