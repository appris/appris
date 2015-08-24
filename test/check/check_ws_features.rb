#!/usr/bin/ruby -w

=begin
  
  Script that download the features files of species from Ensembl FTP

  Send an email when something is wrong
  
=end

# current dir that is included in the search path
BEGIN {
  $LOAD_PATH.unshift(File.dirname(__FILE__))
}

require 'ftpData'
require 'send_email'
require 'yaml'
require 'json'
require 'pp'

if $0 == __FILE__
  error = false
  email_subj = ''
  email_body = ''  
  begin
    # ARGV[0] is the configuration file of species
    config_file = ARGV[0]
        
    spe = File.read(config_file)
    species = JSON.parse(spe)  
  
    # Create ftpData for each species
    species.each do |specie_id,specie|
      ftpdata = FtpData.new(specie, Logger::DEBUG)
      if ftpdata.valid?
        error = false
        email_subj = "The download of files of #{specie_id} were successful"
      else
        error = true
        email_subj = "ERROR downloading ftp files"
        email_body << "TRACE: \n#{ftpdata.strio.string}"
      end      
    end

  rescue Exception => e
      error = true
      email_subj = "ERROR checking the features of species"
      email_body = e.message   
  ensure
    puts "Error: #{error}"
    puts "Subj: #{email_subj}"
    puts "Body: #{email_body}"    
    # Send email
#    if error
#      send_email("jmrodriguez@cnio.es", {
#                                  :subject  => "[APPRIS] #{email_subj}",
#                                  :body     => email_body
#                                })      
#    end
  end
end