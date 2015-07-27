#!/usr/bin/ruby -w

require 'net/smtp'

def send_email(to,opts={})
  opts[:server]      ||= 'webmail.cnio.es'
  opts[:from]        ||= 'jmrodriguez@cnio.es'
  opts[:from_alias]  ||= 'APPRIS Webmaster'
  opts[:subject]     ||= "You need to see this"
  opts[:body]        ||= "Important stuff!"

  msg = <<END_OF_MESSAGE
From: #{opts[:from_alias]} <#{opts[:from]}>
To: <#{to}>
Subject: #{opts[:subject]}

#{opts[:body]}
END_OF_MESSAGE

  Net::SMTP.start(opts[:server]) do |smtp|
    smtp.send_message msg, opts[:from], to
  end
end