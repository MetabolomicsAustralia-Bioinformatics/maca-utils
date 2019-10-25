library("MetaboAnalystR")


call_maca_normalization <- function(fn_auc_csv) {
    "Using MetaboAnalystR, does, in the following order: 
    1. Missing value imputation (replaces all zeros with half of the smallest nonzero value)
    2. Row-wise median normalization
    3. Log2-transformation data
    Implemented because I can't be bothered to re-write median normalization from scratch.


    PARAMS
    ------
    fn_auc_csv: str; path to input csv

    RETURNS
    -------
    tbl0: tibble of normalized data, where first column are sample names and second column are groups. 

    Validated against the metaboanalyst.ca output with sepsis data (10 csvs), with an elementwise error of <10e15. 
    "
    mSet<-InitDataObjects("conc", "stat", FALSE)
    mSet<-Read.TextData(mSet, fn_auc_csv, "rowu", "disc");
    mSet<-SanityCheckData(mSet)
    mSet<-ReplaceMin(mSet);
    mSet<-PreparePrenormData(mSet)
    mSet<-Normalization(mSet, "MedianNorm", "LogNorm", "NULL", ratio=FALSE, ratioNum=20)

    # Load and merge the groups column, which mSet annoyingly doesn't have
    tbl0 <- as_tibble(mSet$dataSet$norm, rownames = "Sample")

    tbl_m <- suppressMessages(read_csv(fn_auc_csv))
    sample_colname <- colnames(tbl_m[,1])
    group_colname <- colnames(tbl_m[,2])
    tbl_m <-  tbl_m %>% dplyr::select(c(!!sample_colname, !!group_colname)) %>% arrange(!!group_colname) %>% dplyr::rename("Sample"=!!sample_colname)
    tbl0 <- inner_join(tbl_m, tbl0, by=!!sample_colname)

    return(tbl0)
}


call_maca_pw_analysis <- function(fn_auc_csv, kegg_species_id) {
	"Calls the pathway enrichment analysis module from MetaboAnalystR. 
	Does row-wise median-normalization and log-transforms the data.
	P-values of pathway enrichment are calculated using the `globaltest` algorithm, and pathway impact scores
	computed using the pathway centrality option. But impact should be disregarded as an overly-abstract
	graph theoretic notion that doesn't necessarily have any biological relevance. 

	PARAMS
	------
	fn_auc_csv: str; filename of input run summary table as a csv file, with AUCs as values. 
	rownames are the sample names, column names are the metabolite names. Column 1 are the 
	experimental groupings. Because of the way this module works, only 2 groups are supported.

	RETURNS
	-------
	list of two outputs:
	tbl.out: output tibble of pathway analysis enrichment. columns:
		metabolite (compound common name), total cmpd, Hits, raw p (raw p value), -log p, 
		Holm adjust(ed p value), FDR, Impact. 
	pw.dict: named list of lists; each key is the pathway ID, and each value is a list of 
		compounds from the input data which appear in that particular pathway. 
	"

	mSet<-InitDataObjects("conc", "pathqea", FALSE)
	mSet<-Read.TextData(mSet, fn_auc_csv, "rowu", "disc");
	mSet<-CrossReferencing(mSet, "name");
	mSet<-CreateMappingResultTable(mSet)
	mSet<-SanityCheckData(mSet)
	mSet<-ReplaceMin(mSet);
	mSet<-PreparePrenormData(mSet)
	mSet<-Normalization(mSet, "MedianNorm", "LogNorm", "NULL", ratio=FALSE, ratioNum=20)
	mSet<-SetKEGG.PathLib(mSet, kegg_species_id)
	mSet<-SetMetabolomeFilter(mSet, F);
	mSet<-CalculateQeaScore(mSet, "rbc", "gt")

	tbl.out <- as_tibble(mSet$analSet$qea.mat, rownames="pw_name")
	pw.dict <- mSet$analSet$qea.hits

	return(list(tbl.out, pw.dict))
}


get_de_metabs <- function(tbl, input_alpha, grp_numerator, grp_denominator, input_fdr) {
    "Takes a tibble of cleaned, transformed values as input, and does:
    1. t-tests at input_alpha
    2. BH correction at input_fdr
    3. Computes FC of averages. Assumes log2 abundances as input, so this is computed as mu1 - mu2. 
    4. Gets FC colours

    PARAMS
    ------
    tbl: input tibble of metabolite abundances, with sample names in col1 and group names in col2. 
    input_alpha: float; alpha at which t-tests are applied.
    grp_numerator: numerator of class/group for fold change calculations
    grp_denominator: denominator of class/group for fold change calculations
    input_fdr: float; threshold p-value at which the benjamini-hochberg method is applied. All results are returned; only affects FC colour assigned. 

    OUTPUT
    ------
    tbl: tibble with the following column names:
        'Sample': str; compound name
        'fc': float; fold change of log2-abundances. 
        'KEGG': str; KEGG Id, retrieved through MetaboAnalystR. 
        'HMDB': str; HMDB Id, retrieved through MetaboAnalystR.
        'raw_p_val': float; raw p-value, output through the t-test.
        'adj_p_val': float; adjusted p-value, computed through the Benjamini-Hochberg procedure.
        'fc_colour': str; colour of fold changes, in hex colour code. 

    NOTES
    -----
     * Unknown KEGG and HMDB Ids get replaced with string `undef`. 
     * Doesn't actually matter if the input tbl has more than 2 groups; groups that are not required do not get
     picked up by %>% filter() anyway. 
     * In FC computations, assumes that the input tibble already has log-abundances, therefore FC = log(x) - log(y), 
     s.t. FC = log(x/y) 
     * Thresholding occurs at >log2(1.25) and <log2(0.75). 
    "
    group.name.col <- names(tbl[,2])
    metab_names <- colnames(tbl)[3:length(colnames(tbl))]
    # Init vec of t.stats p-vals, and FCs
    # Compute t-stats
    p.vals.ls <- vector(mode="numeric", length = length(metab_names))
    fc.ls <- vector(mode="numeric", length = length(metab_names))
    for (i in 1:length(metab_names)) {
        g1.vec <- as.vector(unlist(tbl %>% filter(!!sym(group.name.col)==grp_numerator) %>% select(metab_names[i])))
        g2.vec <- as.vector(unlist(tbl %>% filter(!!sym(group.name.col)==grp_denominator) %>% select(metab_names[i])))

        x <- t.test(g1.vec, g2.vec, conf.level=1-input_alpha)
        p.vals.ls[i] <- x$p.value

        mu1 <- mean(g1.vec)
        mu2 <- mean(g2.vec)
        fc.ls[i] <- mu1 - mu2
    }
    

    names(p.vals.ls) <- metab_names
    names(fc.ls) <- metab_names

    # adjust: BH correction
    adj.p.vals.ls <- p.adjust(p.vals.ls, method = "hochberg", n = length(p.vals.ls))

    # Compute ipath colours: non-significant, FC-up, FC-down, FC-neutral
    fc.colour.ls <- rep("#ACACAC", ncol(tbl)) # default gray (non-significant)
    names(fc.colour.ls) <- metab_names
    for (nm in metab_names) {
        if (adj.p.vals.ls[nm] < input_fdr) {
            if (fc.ls[nm] > 0.22314) {
                fc.colour.ls[nm] <- "#0571b0" #blue
            } else if (fc.ls[nm] < -0.28768) {
                fc.colour.ls[nm] <- "#ca0020" #red
            } else {
                fc.colour.ls[nm] <- "#000000" #black
            }
        }
    }

    # Get all KEGG IDs
    chem.id.tbl <- lookup_chem_id(metab_names)
    kegg.id.vec <- as.vector(unlist(chem.id.tbl %>% select("KEGG")))
    
    # enframe and merge all named lists
    fc.tbl <- tibble::enframe(fc.ls) %>% rename("log2_fc"=value, "Sample"=name)
    p.val.tbl <- tibble::enframe(p.vals.ls) %>% rename("raw_p_val"=value, "Sample"=name)
    adj.p.val.tbl <- tibble::enframe(adj.p.vals.ls) %>% rename("adj_p_val"=value, "Sample"=name)
    fc.colour.tbl <- tibble::enframe(fc.colour.ls) %>% rename("fc_colour"=value, "Sample"=name)
    tbl <- list(fc.tbl, chem.id.tbl, p.val.tbl, adj.p.val.tbl, fc.colour.tbl) %>% reduce(inner_join, by = "Sample")
    
    return(tbl)
}


get_kegg_pw_ref_tbl <- function(kegg_species_id) {
	"Retrieves list of pathways and their associated KEGG Ids. 
	Note that KEGG name mapping can be weird: 
	'Glycolysis / Gluconeogenesis' a.k.a. 'Glycolysis or Gluconeogenesis'
	"
	kl.dt <- as.data.table(keggList("pathway", kegg_species_id), keep.rownames = T)
	kl.dt[, 2] <- as.vector(unlist(kl.dt[ , lapply(V2, function(x) strsplit(x, " - ")[[1]][1])]))
	kl.dt[, 1] <- as.vector(unlist(kl.dt[ , lapply(V1, function(x) gsub("path:", "", x))]))
	colnames(kl.dt) <- c("pw_id", "pw_name")

	return(kl.dt)
}


lookup_chem_id <- function(cpd_names_vec) {
    "Look up various IDs on KEGG, through MetaboAnalystR.
    Keeps only KEGG and HMDB IDs.
    Cleans up output so that blank or <NA> cells become 'undef'
    Requirements: MetaboAnalystR
    
    PARAMS
    ------
    cpd_names_vec: vector of characters; common chemical names of compounds
    
    RETURNS
    -------
    cpd_names_tbl: tibble of IDs.
    "
    # Call to Kegg
    mset <- InitDataObjects("NA", "utils", FALSE)
    mset <- Setup.MapData(mset, cpd_names_vec)
    mset <- CrossReferencing(mset, "name", T, T, T, T, T)
    mset <- CreateMappingResultTable(mset)
    
    # print warnings
    print(mset$msgset$nmcheck.msg[2])
    cpd_names_tbl <- as_tibble(mset$dataSet$map.table)
    cpd_names_tbl <- cpd_names_tbl %>% dplyr::rename("Sample"="Query")
    cpd_names_tbl <- cpd_names_tbl %>% dplyr::select("Sample", "KEGG", "HMDB")

    # replace NA, string "NA", or empty cell with string "undef"
    cpd_names_tbl[cpd_names_tbl == ""] <- "undef"
	cpd_names_tbl[cpd_names_tbl == "NA"] <- "undef"
	cpd_names_tbl <- cpd_names_tbl %>% replace(., is.na(.), "undef")
    
    return(cpd_names_tbl)
}


get_significant_ellipse_coords <- function(url, selection.str, kegg_species_id) {
	"DEPRECATED?
    Makes a POST req to iPath with `selection.str` (which should comprise all
	metabolite hits), and `whole_modules` set to 0. Scrapes the resulting output
	to retrieve the `cx` and `cy` attributes of each `ellipse` html object. 
	This function got modularized for containment. 

	PARAMS
	------
	url: string; API url.
	selection.str: vector of characters, comprising the KEGG Ids of metabs of interest
	kegg_species_id: three letter kegg species id

	RETURNS
	-------
	coords_ls: mtrix of coordinates, shape (n_significant_ellipses, 2)
	"
	selection.str <- paste0(as.vector(unlist(selection.str)), collapse=" W20 #ff0000\n")
	call1 <- list(selection = selection.str, 
              export_type="svg", 
              default_opacity="0.7",
              default_width="1",
              default_radius="5",
              whole_modules="0", 
              tax_filter=kegg_species_id)
	r <- POST(url, body = call1, encode = "form")
	print(http_status(r)$message)
	# Get (cx, cy) coords of significant ellipses, i.e. where fill=<red>
	ellipse_ls <- xml_find_all(read_html(content(r, "text")), '//ellipse')
	x0.attr.dt <- as.data.table(t(as.data.table(xml_attrs(ellipse_ls))))
	colnames(x0.attr.dt) <- names(xml_attrs(ellipse_ls)[[1]])
	x0.attr.dt <- x0.attr.dt[fill=="#ff0000"]
	coords_ls <- c()
	for (i in 1:nrow(x0.attr.dt)) {
	  coords_ls <- rbind(coords_ls, c(x0.attr.dt[i, cx], x0.attr.dt[i, cy]))
	}
	return(coords_ls)
}


clean_and_transform_tibble <- function(tbl, log_bool, z_transform_bool) {
    "DEPRECATION WARNING: Should be deprecated in favour of call_maca_normalization
    "
    cat("DEPRECATION WARNING: Should be deprecated in favour of call_maca_normalization.\n")
    # replace NA with zero
    tbl <- tbl %>% replace(., is.na(.), 0)
    
    # Replace 0 with half of smallest nonzero value
    metab_names <- colnames(tbl)[3:length(colnames(tbl))]
    min.nonzero.val <- sort(unique(as.vector(as.matrix(tbl %>% select(metab_names)))))[2]/2
    tbl[tbl == 0] <- min.nonzero.val
    
    # Take log (ugh), transforming in-place
    if (log_bool) {
        for (nm in metab_names) {
            tbl[nm] <- log2(unlist(tbl[nm]))
        }
    }
    # Do z-transform
    if (z_transform_bool) {
        for (nm in metab_names) {
            tmp <- unlist(tbl[nm])
            tbl[nm] <- (tmp - mean(tmp))/sd(tmp)
        }
    }
    return(tbl)
}


get_ipath_selection_str <- function(tbl, metab_name_colname, kegg_colname, fc_colour_colname, node_width_text) {
    "Get the ipath selection string from an input tibble. Usually goes after get_results_tibble() after the appropriate 
    filter()-ing or select()-ing. 

    PARAMS
    ------
    tbl: input tibble. 
    metab_name_colname: str; name of column containing metabolite names. 
    kegg_colname: str;name of column containing KEGG IDs
    fc_colour_colname: str; name of column containing colours in hex colour code format. 
    node_width_text: str; width of node, for ipath entry. Usually 'W20'

    RETURNS
    -------
    selection_str
    "

    # Extract columns as named lists, named by metab names
    kegg.vec <- as.vector(unlist(tbl[kegg_colname]))
    names(kegg.vec) <- as.vector(unlist(tbl[metab_name_colname]))
    fc.colour.ls <- as.vector(unlist(tbl[fc_colour_colname]))
    names(fc.colour.ls) <- as.vector(unlist(tbl[metab_name_colname]))

    selection_str_ls <- c()
    for (nm in names(kegg.vec)) {
        selection_str_ls <- c(selection_str_ls, paste(kegg.vec[nm], node_width_text, fc.colour.ls[nm], sep=" "))
    }

    selection_str <- paste0(selection_str_ls, collapse = "\n")
    return(selection_str)
}
