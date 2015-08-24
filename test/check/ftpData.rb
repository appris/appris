#!/usr/bin/ruby -w

require 'logger'
require 'stringio'
require 'net/ftp'
require 'fileutils'
require 'zlib'
require 'pp'


# Class that download data files from Ensembl
class FtpData
  
  attr_reader :logger
  attr_reader :strio
  
  def initialize(specie, loglevel=Logger::DEBUG)
    begin
      # logger
      @strio = StringIO.new
      @logger = Logger.new(@strio)
      #$stdout.close
      #$stderr.close
      @logger.level = loglevel
      
      # specie id and specie name
      specie_sci = specie['scientific']
      specie_name = specie_sci.sub! ' ', '_'
      specie_id = specie_sci.downcase
      
      # constant values of remote links (ftp)
      remote = Hash[
        'ftp'  => "ftp.ensembl.org",
        'cdna'  => {
          'dir'  => '/pub/release-__ensembl__/fasta/__specie__id__/cdna',
          'file'  => '__specie__name__.__assembly__.cdna.all.fa.gz'          
        },
        'pep'  => {
          'dir'  => '/pub/release-__ensembl__/fasta/__specie__id__/pep',
          'file'  => '__specie__name__.__assembly__.pep.all.fa.gz'          
        },
        'data'  => {
          'dir'  => '/pub/release-__ensembl__/gtf/__specie__id__',
          'file'  => '__specie__name__.__assembly__.__ensembl__.gtf.gz'          
        }
      ]      
      # constant values of local files
      local = Hash[
        'dir'  => "#{ENV['APPRIS_HOME']}/ws/features/__specie__id__/ensembl__ensembl__",
        'cdna'  => "__specie__id__.transc.fa",
        'pep'   => "__specie__id__.transl.fa",
        'data'  => "__specie__id__.annot.gtf"          
      ]
            
      # replace constant values (id and name)
      self.replace(remote,'__specie__id__', specie_id)
      self.replace(remote,'__specie__name__', specie_name)
      self.replace(local,'__specie__id__', specie_id)
      self.replace(local,'__specie__name__', specie_name)
    
      # replace variable values (assembly and ensembl version)
      specie['assemblies'].each do |ass_id,ensembls|
        assembly,*rest = ass_id.split("|")
        rest = nil
                
        ensembls.each do |ensembl|
          # Clone 
          rem = Marshal.load(Marshal.dump(remote))
          loc = Marshal.load(Marshal.dump(local))
          
          self.replace(rem,'__ensembl__', ensembl.to_s)
          self.replace(rem,'__assembly__', assembly)
          self.replace(loc,'__ensembl__', ensembl.to_s)
          self.replace(loc,'__assembly__', assembly)
          
          unless self.local?(loc) 
            puts ">> SP #{specie_id} > #{assembly} > #{ensembl} does not exit"            
            self.download(rem,loc)
          else
            puts ">> SP #{specie_id} > #{assembly} > #{ensembl} exits"            
          end
        end
      end
      
                  
    rescue Exception => e    
      sms = "#{self.class.name}: "
      sms << "Unable to initialize class: #{e.message}\n"
      sms << e.backtrace.inspect
      @logger.error(sms)
    end
  end
  
  # everything was ok
  def valid?()
    val = nil
    begin
      if @strio.size == 0
        val = true
      else
        val = false
      end
    end
    return val
  end  
  
protected

  # Replace values
  def replace(data,old,new)
    begin
      data.each do |key,val|
        if val.is_a?(Hash)
          val.each do |k,v|
            v.gsub! old, new
          end
        else
          val.gsub! old, new
        end
      end      
    rescue Exception => e
      sms = "#{self.class.name}: "
      sms << e.backtrace.inspect
      @logger.error(sms)    
    end
  end
  
  # Have we got already the species files  
  def local?(loc)
    local = nil
    begin
      loc.each do |key,val|
        if ( key == 'cdna' || key == 'pep' || key == 'data' )
          f = "#{loc['dir']}/#{val}"
          if ( !File.file?(f) )
            local = false
            return local
          else
            local = true
          end
        end
      end
    rescue Exception => e
      sms = "#{self.class.name}: "
      sms << e.backtrace.inspect
      @logger.error(sms)    
    end    
    return local
  end  
  
  # Download files using ftp
  def download(rem,loc)
    begin
      
      # create local dir      
      FileUtils.mkdir_p loc['dir']
        
      # ftp connection        
      Net::FTP.open(rem['ftp']) do |ftp|
        ftp.login('anonymous','')        
        rem.each do |key,r_val|
          if ( key == 'cdna' || key == 'pep' || key == 'data' )
            r_dir = r_val['dir']
            r_file = r_val['file']
            l_file = "#{loc['dir']}/#{loc[key]}.gz"
            @logger.info(" #{r_dir}/#{r_file} > #{l_file}")
            ftp.chdir(r_dir)
            ftp.getbinaryfile(r_file,l_file)
            self.uncompress(l_file)
          end
        end 
      end
              
    rescue Exception => e
      sms = "#{self.class.name}: "
      sms << e.backtrace.inspect
      @logger.error(sms)    
    end    
  end
   
  def uncompress(file)
    begin
      #Archive::Tar.extract(file)
      #gz = Zlib::GzipReader.new(file) 
      #result = gz.read
      #puts result      
      system("gzip -d #{file}")
      rescue Exception => e
        sms = "#{self.class.name}: "
        sms << e.backtrace.inspect
        @logger.error(sms)    
    end    
  end
    
end
