ALTER TABLE CONSENSUS ADD INDEX (`PDBID`);
ALTER TABLE CONSENSUS ADD INDEX (`CLUSTID`);

ALTER TABLE INFOACC ADD INDEX (`EC1`);
ALTER TABLE INFOACC ADD INDEX (`EC2`);
ALTER TABLE INFOACC ADD INDEX (`EC3`);
ALTER TABLE INFOACC ADD INDEX (`EC4`);
ALTER TABLE INFOACC ADD INDEX (`PDBID`);
ALTER TABLE INFOACC ADD INDEX (`UNIACC1`);
ALTER TABLE INFOACC ADD INDEX (`UNIACC2`);
ALTER TABLE INFOACC ADD INDEX (`UNIACC3`);

ALTER TABLE COMPARE35 ADD INDEX (`CLUSTID`);
ALTER TABLE COMPARE40 ADD INDEX (`CLUSTID`);
ALTER TABLE COMPARE45 ADD INDEX (`CLUSTID`);
ALTER TABLE COMPARE35 ADD INDEX (`TEMPID`);
ALTER TABLE COMPARE40 ADD INDEX (`TEMPID`);
ALTER TABLE COMPARE45 ADD INDEX (`TEMPID`);
ALTER TABLE COMPARE35 ADD INDEX (`CLUSTSITEID`);
ALTER TABLE COMPARE40 ADD INDEX (`CLUSTSITEID`);
ALTER TABLE COMPARE45 ADD INDEX (`CLUSTSITEID`);
ALTER TABLE COMPARE35 ADD INDEX (`TEMPSITEID`);
ALTER TABLE COMPARE40 ADD INDEX (`TEMPSITEID`);
ALTER TABLE COMPARE45 ADD INDEX (`TEMPSITEID`);

ALTER TABLE BINDSITE35 ADD INDEX (`PDBID`);
ALTER TABLE BINDSITE40 ADD INDEX (`PDBID`);
ALTER TABLE BINDSITE45 ADD INDEX (`PDBID`);
ALTER TABLE BINDSITE35 ADD INDEX (`CLUSTID`);
ALTER TABLE BINDSITE40 ADD INDEX (`CLUSTID`);
ALTER TABLE BINDSITE45 ADD INDEX (`CLUSTID`);
ALTER TABLE BINDSITE35 ADD INDEX (`CADID`);
ALTER TABLE BINDSITE40 ADD INDEX (`CADID`);
ALTER TABLE BINDSITE45 ADD INDEX (`CADID`);

ALTER TABLE CSITE35 ADD INDEX (`CLUSTID`);
ALTER TABLE CSITE40 ADD INDEX (`CLUSTID`);
ALTER TABLE CSITE45 ADD INDEX (`CLUSTID`);
ALTER TABLE CSITE35 ADD INDEX (`PDBID`);
ALTER TABLE CSITE40 ADD INDEX (`PDBID`);
ALTER TABLE CSITE45 ADD INDEX (`PDBID`);

ALTER TABLE SITE35 ADD INDEX (`BINDID`);
ALTER TABLE SITE40 ADD INDEX (`BINDID`);
ALTER TABLE SITE45 ADD INDEX (`BINDID`);
ALTER TABLE SITE35 ADD INDEX (`CADID`);
ALTER TABLE SITE40 ADD INDEX (`CADID`);
ALTER TABLE SITE45 ADD INDEX (`CADID`);

ALTER TABLE CCTEVAL_35 ADD INDEX (`CLUSTID`);
ALTER TABLE CCTEVAL_40 ADD INDEX (`CLUSTID`);
ALTER TABLE CCTEVAL_45 ADD INDEX (`CLUSTID`);

