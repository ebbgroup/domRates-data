/*
 * DomainAll.hpp
 *
 *  Created on: Aug 18, 2015
 *      Author: ckeme_01
 */

#ifndef SRC_DOMAINMODULE_HPP_
#define SRC_DOMAINMODULE_HPP_




/** @defgroup DomainGroup Domain module
 * \brief Module containing all domain related classes and functions.
 *
 *  Module with all classes & functions related to domains, including classes for domains, domain arrangements and domain arrangement sets. Three types of domains are currently implemented:
 *  Domain, the base class and two derived classes (PfamDomain and SFDomain) implementing special features of the Pfam and the Superfamily databases. A DomainArrangement represents the
 *  order of domains found in a protein. The DomainArrangementSet is a collection of domain arrangements, e.g. all domain arrangements of a proteome.
 *
 *	Below you can find some example for reading and accessing domains.
 *  \code
 *  // reading a DomainArrangementSet
 *  BioSeqDataLib::DomainArrangementSet<BioSeqDataLib::PfamDomain> arrangementSet;
 *	arrangementSet.read("path/to/pfamscan.pl/output");
 *
 *	// accessing the domain arrangement of a specific sequence
 *	BioSeqDataLib::DomainArrangement<BioSeqDataLib::PfamDomain> da = arrangmentSet["sequence name"]
 *
 *	// accessing a specific domain in the arrangement
 *	BioSeqDataLib::PfamDomain domain = da[<domain position>];
 *
 * \endcode
 *
 *  @{
 */



/**
 * \file DomainModule.hpp
 * \brief Simple header file including all Domain related headers.
 *
 */


#include "domain/Domain.hpp"
#include "domain/DomainExt.hpp"
#include "domain/PfamDomain.hpp"
#include "domain/SFDomain.hpp"
#include "domain/DomainArrangement.hpp"
#include "domain/DomainArrangementSet.hpp"


/** @} */ // Domain module


#endif /* SRC_DOMAINMODULE_HPP_ */
