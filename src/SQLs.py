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
    
    'test': """
    SELECT *
    FROM ComputedFeatures
    WHERE ComputedFeatures.ek_gene_id = ?
    """
}
