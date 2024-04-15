sqlite_object = 'file:.dbcache/EnrichKitDB.sqlite?mode=ro'

SQLs = {
    'ekid2geneid': """
    SELECT *
    FROM ID_Mapper
    WHERE ID_Mapper.ek_gene_id = ?
    """,
    
    'posi2gene': """
    SELECT ek_gene_id
    FROM Gene 
    JOIN Species ON Gene.species = Species.ek_species 
    WHERE Species.name_short = ? AND Gene.seqname = ? AND Gene.start <= ? AND Gene.end >= ?
    """,
    
    'posi2exon': """
    SELECT exon.*
    FROM Exon
    WHERE Exon.ek_gene_id = ?
    AND Exon.start <= ?
    AND Exon.end >= ?
    """,
    
    'getGenelimit': """
    SELECT Gene.start,Gene.end, Gene.strand, Gene.gene_biotype
    FROM Gene
    WHERE Gene.ek_gene_id = ?
    """,
    
    'posi2feature': """
    SELECT *
    FROM Feature
    WHERE Feature.ek_gene_id = ?
    AND Feature.start <= ?
    AND Feature.end >= ?
    """,
    
    'posi2cfeature1': """
    SELECT *
    FROM ComputedFeatures
    WHERE ComputedFeatures.ek_gene_id = ?
    AND ComputedFeatures.start <= ?
    AND ComputedFeatures.end >= ?
    """,
    
    'posi2cfeature_sd': """
    SELECT *
    FROM ComputedFeatures
    WHERE ComputedFeatures.ek_gene_id = ?
    AND ComputedFeatures.feature == 'splice donor'
    AND ComputedFeatures.start <= ?
    AND ComputedFeatures.end >= ?
    """,
    
    'posi2cfeature_sa': """
    SELECT *
    FROM ComputedFeatures
    WHERE ComputedFeatures.ek_gene_id = ?
    AND ComputedFeatures.feature == 'splice acceptor'
    AND ComputedFeatures.start <= ?
    AND ComputedFeatures.end >= ?
    """,

    'id_convert_ensembl': """
    SELECT ID_Mapper.*
    FROM ID_Mapper
    JOIN Species ON ID_Mapper.species = Species.ek_species
    WHERE Species.name_short = ? AND ID_Mapper.gene_id = ?
    """,

    'id_convert_entrez': """
    SELECT ID_Mapper.*
    FROM ID_Mapper
    JOIN Species ON ID_Mapper.species = Species.ek_species
    WHERE Species.name_short = ? AND ID_Mapper.entrez_id = ?
    """,

    'extract_geneset_ensembl': """
    SELECT FilteredPathway.pathway_id, FilteredPathway.pathway_description, ID_Mapper.gene_id
    FROM Involve
    JOIN ID_Mapper ON Involve.ek_gene_id = ID_Mapper.ek_gene_id
    JOIN (
        SELECT Pathway.ek_pathway_id, Pathway.pathway_id, Pathway.pathway_description
        FROM Pathway
        JOIN Pathway_Meta ON Pathway.pathway_meta = Pathway_Meta.pathway_meta_id
        WHERE Pathway_Meta.species = ? AND Pathway_Meta.name = ?
    ) AS FilteredPathway ON Involve.ek_pathway_id = FilteredPathway.ek_pathway_id;
    """,

    'extract_geneset_entrez': """
    SELECT FilteredPathway.pathway_id, FilteredPathway.pathway_description, ID_Mapper.entrez_id,
    FROM Involve
    JOIN ID_Mapper ON Involve.ek_gene_id = ID_Mapper.ek_gene_id
    JOIN (
        SELECT Pathway.ek_pathway_id, Pathway.pathway_id, Pathway.pathway_description
        FROM Pathway
        JOIN Pathway_Meta ON Pathway.pathway_meta = Pathway_Meta.pathway_meta_id
        WHERE Pathway_Meta.species = ? AND Pathway_Meta.name = ?
    ) AS FilteredPathway ON Involve.ek_pathway_id = FilteredPathway.ek_pathway_id;
    """,

    'extract_tf_gene':"""
    SELECT *
    FROM TF_Gene;
    """,

    'test': """
    SELECT *
    FROM ComputedFeatures
    WHERE ComputedFeatures.ek_gene_id = ?
    """
}
