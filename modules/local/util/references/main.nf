def get_genome_attribute(params, attribute) {
    if (params.genomes && params.host_genome && params.genomes.containsKey(params.host_genome)) {
        if (params.genomes[ params.host_genome ].containsKey(attribute)) {
            return params.genomes[ params.host_genome ][ attribute ]
        }
    }
    return null
}
