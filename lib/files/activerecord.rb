
module GenomeDB	
	require 'bio'
	
	class Version < DBConnection
	
	end
		
	class Dna < DBConnection
		has_many :genes
    	has_many :transcripts, :through => :genes
    
    	def sequence
    		path = GenomeDB::Version.find(:first).file_path
    		return Bio::FastaFormat.open("#{path}/fasta/#{self.connection.current_database}.fasta").find{|e| e.definition == self.accession }.naseq
    	end
      
	end

	class Gene < DBConnection
		set_primary_key 'id'
   	 	belongs_to :dna, :foreign_key => "dna_id"
    	has_many :transcripts
    	has_many :xref_gene_external_dbs
    	has_many :xref_gene_samples
      
    	def slice
    		return Slice.new(self.dna,self.start,self.stop,self.strand)
    	end
      
   		def seq
    		return self.slice.seq
    	end
      
      	def longest_transcript        
        	transcripts = self.transcripts
         	if transcripts.length == 1
            	return transcripts.shift
         	else
            answer = nil
            self.transcripts.each do |t|
               if answer.nil?
                  answer = t 
               else
                  answer = t if answer.cds.length < t.cds.length
               end
           	end
         	return answer
         end   
            
      end
      
    end

    class Transcript < DBConnection
  		set_primary_key 'id'
		belongs_to :gene, :foreign_key => "gene_id"
      	has_one :dna, :through => :gene
      	has_many :xref_transcript_exons
      	has_many :exons, :through => :xref_transcript_exons, :order => "start ASC"
      	has_one :translation
      
      	def slice
        	return Slice.new(self.gene.dna,self.start,self.stop,self.strand)
      	end

		def exons
			return self.xref_transcript_exons.collect{|x| x.exon }.sort_by{|e| e.start }
		end
		
      	def seq
        	return self.slice.seq
      	end
      
      	def translation_product
             
        	if self.strand == 1
        		return Bio::Sequence::NA.new(self.cds).translate(self.translation.codon_start,self.translation.transl_table)
        	else
          		return Bio::Sequence::NA.new(self.cds).complement.translate(self.translation.codon_start,self.translation.transl_table)
        	end
       
      	end
      
      	def cds
      		seq = self.gene.dna.sequence
        	cds = ""
        	coding_start,coding_end = self.translation.coding_start,self.translation.coding_end
        	self.exons.each do |exon|
            	if coding_start > exon.start and coding_start < exon.stop and coding_end > exon.start and coding_end < exon.stop
               		cds += seq.subseq(coding_start,coding_end)
            	elsif coding_start > exon.start and coding_start < exon.stop
               		cds += seq.subseq(coding_start,exon.stop)
            	elsif coding_end > exon.start and coding_end < exon.stop
               		cds += seq.subseq(exon.start,coding_end)
            	else
               		cds += seq.subseq(exon.start,exon.stop)
            	end
        	end
        	return cds
      	end
    
	end

    class XrefTranscriptExon < DBConnection
    	set_primary_keys :transcript_id, :exon_id
      	belongs_to :transcript #, :foreign_key => "transcript_id"
      	belongs_to :exon #, :foreign_key => "exon_id"
    end

    class Exon < DBConnection
    	set_primary_key 'id'
    	has_many :xref_transcript_exons
      	has_many :transcripts, :through => :xref_transcript_exons
    end

    class Translation < DBConnection
    	belongs_to :transcript, :foreign_key => "transcript_id"
      
      	def seq
        	self.transcript.translation_product
      	end
      
    end 
    
    class ExternalDb < DBConnection
    	has_many :xref_gene_external_dbs
    end
    
    class XrefGeneExternalDb < DBConnection
    	belongs_to :gene, :foreign_key => "gene_id"
    	belongs_to :external_db, :foreign_key => "external_db_id"
    end
    
    class Sample < DBConnection
    	has_many :xref_gene_samples
    end
    
    class XrefGeneSample < DBConnection
    	belongs_to :gene, :foreign_key => "gene_id"
    	belongs_to :sample, :foreign_key => "sample_id"
    end
    	
    class Slice
      
      	attr_reader :dna, :start, :stop, :strand
      	def initialize(dna,start,stop,strand)
        	@dna,@start,@stop,@strand = dna,start,stop,strand
      	end
      
      	def seq
        	answer = nil
        	self.strand == 1 ? answer = self.dna.sequence.subseq(start,stop) : answer = self.dna.sequence.subseq(start,stop).complement
        	return answer
      	end
      
    end

end # GenomeDB
