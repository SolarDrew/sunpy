# -*- coding: utf-8 -*-

import os
import numpy as np
import datetime as dt
from glob import glob


def aia_bp_read_response_table(tablefile, silent=False):
    response_dir = '~sunpy/data/'
    if len(tablefile) == 0:
        tablefile = os.join(response_dir, 'aia_V2_response_table.txt')

    hastable = os.path.isfile(tablefile)
    if hastable:
        #######temp_tablefile = STRJOIN(STRSPLIT(/extract, tablefile, '*'))
        temp_tablefile = tablefile.replace('*', '')
        tabledat = RD_TFILE(temp_tablefile, 14)
    else:
        if not silent:
            print 'aia_bp_read_response_table: Missing file ' + tablefile
        return 0

    table_template = {'date' : '', 't_start': '', 't_stop': '', 'ver_num': 0,
                      'wave_str': '', 'wavelnth': 0, 'eperdn': 0.0,
                      'dnperpht': 0.0, 'eff_area': 0.0, 'eff_wvln': 0.0,
                      'effa_p1': 0.0, 'effa_p2': 0.0, 'effa_p3': 0.0,
                      'rmse' : 0.0}

    tablesize = tabledat.size()
    for i in range(1, tablesize[2]):
        thisline = table_template.copy()
        thisline['date'] = tabledat[0,i]
        thisline['t_start'] = tabledat[1,i]
        thisline['t_stop'] = tabledat[2,i]
        thisline['ver_num'] = tabledat[3,i]
        thisline['wave_str'] = tabledat[4,i]
        thisline['wavelnth'] = tabledat[5,i]
        thisline['eperdn'] = tabledat[6,i]
        thisline['dnperpht'] = tabledat[7,i]
        thisline['eff_area'] = tabledat[8,i]
        thisline['eff_wvln'] = tabledat[9,i]
        thisline['effa_p1'] = tabledat[10,i]
        thisline['effa_p2'] = tabledat[11,i]
        thisline['effa_p3'] = tabledat[12,i]
        thisline['rmse'] = tabledat[13,i]
    if i == 1:
        tablestr = thisline
    else:
        tablestr = [tablestr, thisline]

    return tablestr


def aia_bp_parse_effarea(oldfullresp, sampleutc=None, resptable=None, 
                         version=None, evenorm=None, timedepend=None,
                         ver_date=None, channels=None, dn=False):#, uv=uv):
    if len(version) == 0:
        version = 2
    otags = np.array(oldfullresp.keys())
    #ntags = len(oldfullresp)

    #+++++++++++++++++++++++++++++++++
    # Add the informational fields from the old response structure to the new 
    # response structure
    #--------------------------------
    newfullresp = {'name': oldfullresp.name,
                   'version': int(version), # FIX?
                   'date': sampleutc, #AIA_BP_UTC2DATE_STRING(sampleutc), # WHAT?!
                   'filedate': oldfullresp['date']}

    #+++++++++++++++++++++++++++++++++
    # Decide which channels will be included, and set up the field
    # appropriately
    #--------------------------------
    if channels:
        cindex = np.where_ARR(oldfullresp['channels'], channels.upper()) # ??
    else:
        channels = oldfullresp['channels']
        cindex = range(len(channels)) # Can't remember what indgen does
    newfullresp['channels'] = oldfullresp['channels'][cindex]

    #++++++++++++++++++++++++++++++++
    # Generate substructure describing corrections to the instrument 
    # calibration
    #--------------------------------
    if len(resptable) > 0:
        corr_str = AIA_BP_CORRECTIONS(resptable, sampleutc, oldfullresp, ver_date, effarea=True, uv=uv) # ??? AAAAAAAARGH
        if evenorm:
            corr_str.evenorm_applied = 'YES' # ???
            evechan = corr_str.evenorm['channels'][:3]
            evescale = corr_str.evenorm.aia_over_eve # ???
        if timedepend:
            corr_str.time_applied = 'YES' # ??
            timestr = corr_str.time # ??
            timechan = timestr['channels'][:3]
            timescale = timestr.ea_sample_over_t0 # ??
    else:
        # For some reason, the response table information was not passed in. 
        # Silently ignore all requests to include corrections
        evechan = np.array(channels)
        evescale = np.array(len(channels), type=np.float64) + 1.0
        timechan = np.array(channels)
        timescale = np.array(len(channels), type=np.float64) + 1.0
        if evenorm or timedepend:
            print 'aia_bp_parse_effarea: No input response table', 'Cannot apply evenorm or timedepend corrections'
            evenorm = 0
            timedepend = 0

    #+++++++++++++++++++++++++++++++++
    # Loop through each channel and add it to the new response structure, 
    # correcting the units and format of the substructure as you go
    #--------------------------------
    for i in range(len(cindex)-1):
        thischan = channels[i].upper()
        thischanstr = None #oldfullresp.(np.where(otags == thischan))
        thisfullstr = None #STR_SUBSET(oldfullresp.(np.where(otags == thischan)+1), 'name,date,wave,effarea,units,geoarea,scale,platescale,wavemin,wavestep,wavenumsteps,filtersincludemesh,contamthick,usecontam,usephottoelec,usephottodn,elecperev,elecperdn,useerror,extrapolatecomponents,ent_filter,fp_filter,primary,secondary,ccd,contam,cross_area')
        if i == 0:
            wave = thischanstr['wave']
            platescale = thischanstr['platescale']
            units = thischanstr['units']
            smallresp = {'name': newfullresp.name,
                         'version': newfullresp.version,
                         'date': newfullresp.date, 'channels': channels,
                         'wave': wave}
        effarea = thisfullstr['effarea']
        if dn:
            evperphot = 12398. / wave
            elecperphot = (evperphot * thisfullstr.elecperev) > 1
            elecperdn = thisfullstr.elecperdn
            effarea = effarea * elecperphot / elecperdn
            units = 'cm^2 DN phot^-1'
            thisfullstr.usephottoelec = 1
            thisfullstr.usephottodn = 1
        if evenorm:
            evechan_id = np.np.where(evechan == thischan[:3])
            this_evescale = evescale[evechan_id]
            effarea = effarea * this_evescale[0]
        if timedepend:
            timechan_id = np.np.where(timechan == thischan[:3])
            this_timescale = timescale[timechan_id]
            effarea = effarea * this_timescale[0]
        shortstr = {'name': thischanstr['name'], 'wave': wave, 'ea': effarea,
                    'platescale': platescale, 'units': units}
        thisfullstr['effarea'] = effarea
        thisfullstr['units'] = units

        if i == 0:
            all = effarea 
            smallrespstr = {thischan, shortstr}
        else:
            all = [[all], [effarea]]
            smallrespstr[thischan] = shortstr

        newfullresp[thischan] = shortstr
        newfullresp[thischan+'_full'] = thisfullstr

    smallresp['all'] = TRANSPOSE(all)
    smallresp['platescale'] = platescale
    smallresp['units'] = units
    #smallrespstr) # ?????

    # Add notes and corrections to the end of the structures
    if len(resptable) > 0:
        newfullresp['corrections'] = corr_str
        smallresp['corrections'] = corr_str
    np.wherenotes = np.np.where(otags == 'NOTES')#, hasnotes)
    if len(np.wherenotes) > 0:#hasnotes > 0:
        newfullresp['notes'] = oldfullresp['notes']

    return newfullresp


def aia_bp_blend_channels(infullresp, short=False):
    """Takes in a long-form instrument response structure and adjusts the primary
    and secondary mirror reflectivity (and, thus, effective area) for the channels
    on telescope 1 and telescope 4 to account for potential cross-contamination
    between the channels.
    """
    resptags = infullresp.keys()
    #wvls = int(infullresp['channels'][1:3])
    wvls = np.array([int(chan[1:3] for chan in infullresp['channels'])])
    w94 = np.where(wvls == 94)#, num94)
    num94 = len(w94)
    w131 = np.where(wvls == 131)#, num131)
    num131 = len(w131)
    w304 = np.where(wvls == 304)#, num304)
    num304 = len(w304)
    w335 = np.where(wvls == 335)#, num335)
    num335 = len(w335)

    outfullresp = infullresp
    if num94 > 0:
        for i in range(num94):
            if num304 < 0:
                return #STOP # Can't calculate blend if you don't have data on the other channel
            thischan = infullresp.channels[w94[i]]
            thisfullindex = np.np.where(resptags == thischan + '_FULL')
            thisindex = np.np.where(resptags == thischan)
            thisstr = infullresp[thisfullindex]
            crosschan = infullresp.a304_full
            crosstalk_effarea = crosschan.ent_filter * crosschan.primary * crosschan.secondary * crosschan.contam * crosschan.ccd * crosschan.geoarea * thisstr.fp_filter
            effarea = thisstr.effarea + crosstalk_effarea
            outfullresp[thisfullindex]['cross_area'] = crosstalk_effarea
            outfullresp[thisfullindex][effarea] = effarea
            outfullresp[thisindex]['ea'] = effarea
    if num131 > 0:
        for i in range(num131):
            if num335 < 0:
                return #STOP	# Can't calculate blend if you don't have data on the other channel
            thischan = infullresp.channels[w131[i]]
            thisfullindex = np.np.where(resptags == thischan + '_FULL')
            thisindex = np.where(resptags == thischan)
            thisstr = infullresp[thisfullindex]
            crosschan = infullresp['a335_full']
            crosstalk_effarea = crosschan['ent_filter'] * crosschan['primary']\
                * crosschan['secondary'] * crosschan['contam'] \
                * crosschan['ccd'] * crosschan['geoarea'] * thisstr['fp_filter']
            effarea = thisstr['effarea'] + crosstalk_effarea
            outfullresp[thisfullindex]['cross_area'] = crosstalk_effarea
            outfullresp[thisfullindex][effarea] = effarea
            outfullresp[thisindex]['ea'] = effarea
    if num304 > 0:
        for i in range(num304):
            if num94 < 0:
                return #STOP	# Can't calculate blend if you don't have data on the other channel
            thischan = infullresp['channels'][w304[i]]
            thisfullindex = np.where(resptags == thischan + '_FULL')
            thisindex = np.where(resptags == thischan)
            thisstr = infullresp[thisfullindex]
            crosschan = infullresp['a94_full']
            crosstalk_effarea = crosschan['ent_filter'] * crosschan['primary']\
                * crosschan['secondary'] * crosschan['contam'] \
                * crosschan['ccd'] * crosschan['geoarea'] * thisstr['fp_filter']
            effarea = thisstr[effarea] + crosstalk_effarea
            outfullresp[thisfullindex]['cross_area'] = crosstalk_effarea
            outfullresp[thisfullindex][effarea] = effarea
            outfullresp[thisindex]['ea'] = effarea
    if num335 > 0:
        for i in range(num335):
            if num131 < 0:
                return #STOP	# Can't calculate blend if you don't have data on the other channel
            thischan = infullresp['channels'][w335[i]]
            thisfullindex = np.where(resptags == thischan + '_FULL')
            thisindex = np.where(resptags == thischan)
            thisstr = infullresp[thisfullindex]
            crosschan = infullresp['a131_full']
            crosstalk_effarea = crosschan['ent_filter'] * crosschan['primary']\
                * crosschan['secondary'] * crosschan['contam'] \
                * crosschan['ccd'] * crosschan['geoarea'] * thisstr['fp_filter']
            effarea = thisstr.effarea + crosstalk_effarea
            outfullresp[thisfullindex]['cross_area'] = crosstalk_effarea
            outfullresp[thisfullindex][effarea] = effarea
            outfullresp[thisindex]['ea'] = effarea

    outfullresp['notes'] = outfullresp['notes'] + ' BLENDED'
    """if short:
        presp = aia_bp_parse_effarea(outfullresp, retval)
    else:
        retval = outfullresp"""

    return outfullresp#retval


def aia_bp_parse_tresp(oldfullresp, smallresp, chiantifix=None, channels=None):
    """Helper routine called by AIA_GET_RESPONSE to reformat temperature response
    structure, select a subset of channels, etc.
    INPUTS:
        oldfullresp	--	Full response structure as read in from .genx file
    OUTPUTS:
    	smallresp	--	Reduced response structure 
    RETURNS:
        newfullresp	--	Full response structure with specified channels, units, reduced
    				set of keywords, etc.
    KEYWORDS:
        channels	--	string array with the names of the channels to be included in
    				the output response structures
        chiantifix	--	[OPTIONAL] a structure containing the T response for channels np.where the
    				CHIANTI data is incomplete or erroneous
    """
    otags = oldfullresp.keys()
    #ntags = len(oldfullresp)

    if len(channels) > 0:
        cindex = np.where_ARR(oldfullresp.channels, channels.upper())
    else:
        channels = oldfullresp.channels
        cindex = range(len(channels))

    newfullresp = STR_SUBSET(oldfullresp, 'name,date,effarea_version,emiss_version,emissinfo,logte,wave,tunits,twunits')
    newfullresp['channels'] = channels.upper()
    newfullresp['tresp'] = oldfullresp['tresp'][:,cindex]
    newfullresp['twresp'] = oldfullresp['twresp'][:,:,cindex]

    smallresp = STR_SUBSET(newfullresp, 'name,date,effarea_version,emiss_version,emissinfo,channels')
    smallresp['units'] = newfullresp['tunits']
    smallresp['logte'] = newfullresp['logte']

    for i in range(len(channels)):
        thischan = channels[i].upper()
        chanindex = np.where(otags == thischan)
        thisstr = oldfullresp[chanindex]
        newfullresp[thischan] = thisstr
        smallrespstr = STR_SUBSET(thisstr, 'name,units,logte,tresp')
        if i == 0:
            all = thisstr['tresp']
            #smallresps[thischan] = smallrespstr
        else:
            all = [[all], [thisstr.tresp]]
            #smallresps[thischan] = smallrespstr

    #+++++++++++++++++++++++++++++++
    # Apply (or don't) the chianti fix (adjusting the 94 and 131 channel temperature response
    # to account for missing lines)
    #-------------------------------
    corrtag = np.where(otags == 'CORRECTIONS')#, hastag)
    if len(corrtag) > 0:
        tgrid = newfullresp['logte']
        corr_str = oldfullresp['corrections']
        if len(chiantifix) > 0:
            chianti_str = chiantifix
            ch_channels = chianti_str.channels
            # If necessary, apply time correction to chianti fixes
            if corr_str.time_applied > 'YES':
                for i in range(len(channels)):
                    chianti_fix_id = np.where(ch_channels == channels[i])#, has_chfix)
                    time_fix_id = np.where(corr_str.time.channels == channels[i])#, has_tfix)
                    if (len(chianti_fix_id) == 1) and (len(time_fix_id) == 1):
                        chianti_str.empirical_minus_raw[:,i] = chianti_str.empirical_minus_raw[:,chianti_fix_id] * corr_str.time.ea_sample_over_t0[time_fix_id]
            # Apply chianti fixes to the temperature response
            nftags = newfullresp.keys()
            #nstags = smallresps.keys()
            for i in range(len(channels)):
                thischan = channels[i].upper()
                fchanindex = np.where(nftags == thischan)
                #schanindex = np.where(nstags == thischan)
                thisstr = newfullresp[fchanindex]
                this_orig_tresp = thisstr.tresp
                chianti_fix_id = np.where(ch_channels == thischan)#, has_chfix)
                if chianti_fix_id.size() == 0:#(not has_chfix):
                    this_fix_tresp = this_orig_tresp
                else:
                    this_fix_tresp = this_orig_tresp + chianti_str.empirical_minus_raw[:,chianti_fix_id]
                newfullresp[fchanindex]['tresp'] = this_fix_tresp
                all[:,i] = this_fix_tresp
                #smallresps[schanindex]['tresp'] = this_fix_tresp
            chfix_applied = 'YES'
        else:
		chfix_applied = 'NO'
		chianti_str = {'EMPIRICAL_MINUS_RAW': np.zeros((len(tgrid), len(channels)), dtype=np.float64),
                          'LOGTE': tgrid, 'CHANNELS': channels}
        ainfo4 = 'CHIANTIFIX: Sub-structure giving a description of the empirical corrections to the CHIANTI emissivity'
        ainfo5 = 'CHIANTIFIX_APPLIED: Does the temperature response use the empirical correction to CHIANTI emissivity?'
        corr_str['CHIANTIFIX'] = chianti_str
        corr_str['CHIANTIFIX_APPLIED'] = chfix_applied
        corr_str['info'][4] = ainfo4
        corr_str['info'][5] = ainfo5
        
        # Add corrections structure
        #smallresp = CREATE_STRUCT(smallresp, 'all', all, smallresps, 'corrections', corr_str)
        newfullresp ['corrections'] = corr_str
    #else:
    #    # For some reason (e.g. you are looking at old data) there is no corrections
    #    # substructure, so just go on without it
    #    smallresp = CREATE_STRUCT(smallresp, 'all', all, smallresps)
    return newfullresp


def aia_get_response(effective_area=False, area=False, temp=False, emiss=False,
	full=False, all=False, uv=False, dn=False, phot=False, noblend=False,
	evenorm=False, timedepend=False, chiantifix=False, version=4,
     emversion=False, respversion=False, silent=False, loud=False):

    aiaresp = '~/sunpy/data/'
    if aiaresp == '':
        aiaresp = os.join('$SSW_AIA','response') # Don't know what $SSW_AIA should be replaced with. Also, there is no response folder
    #aiaresp = (FILE_SEARCH(aiaresp))[0] # Don't know what this is

    area = area or effective_area
    if len(evenorm) > 0:
        user_evenorm = 1-np.array(evenorm)
    else:
        user_evenorm = 0 #	Did the user explicitly set it to evenorm=0?

    #++++++++++++++++++++++++++
    # Detect and warn of incompatible keyword settings
    #--------------------------
    if (temp and emiss and area):
        print 'Please specify only one of the following:'
        print '  temp  emiss  area'
        print 'Defaulting to area...'
        temp = False
        emiss = False
        area = True
    if (not temp and not emiss and not area):
        area = True # Silently switch on /area
    if uv:
        area = False # Henceforth, (area=1) implies EUV area
        all = True # Always return all UV channels
    if temp and uv:
        print 'Cannot generate temperature response for UV channels'
        print 'Returning UV channel effective area'
        temp = False
    if emiss and uv:
        print 'No emissivity defined for UV channels'
        print 'Returning UV channel effective area'
        emiss = False
    if emiss and (dn or noblend or phot or timedepend or evenorm):
        print 'dn, phot, time, evenorm and noblend keywords are not defined for emissivity structure'
        dn = False
        noblend = False
        phot = False
        timedepend = False
        evenorm = False
    if (not (dn and phot)) and (not emiss): # Units not explicitly specified - check with user
        print 'AIA_GET_RESPONSE...'
        print 'Do you want to include the DN/phot factor in the response function?'
        print 'Specify either dn=True or phot=True to avoid this message'
        reply = raw_input('Include DN/phot correction? (Default is yes) [y/n] ')
        if reply.lower() in ['n', 'no', 'nein']:
            print 'AIA_GET_RESPONSE: NOT including DN/phot factor'
            dn = False
            phot = True
        else:
            print 'AIA_GET_RESPONSE: including DN/phot factor'
            dn = True
            phot = False
    if timedepend and not evenorm and not user_evenorm: # Timedepend but not evenorm - check if user is sure
        print 'AIA_GET_RESPONSE...'
        print 'You asked for time dependent correction but not EVE normalization'
        print 'Are you sure you do not want to include EVE normalization?'
        reply = raw_input('Include EVE normalization as well? (Default is no) [y/n] ')
        if reply.lower() in ['y', 'yes']:
            print 'AIA_GET_RESPONSE: Including EVE normalization'
            evenorm = True
        else:
            print 'AIA_GET_RESPONSE: Not including EVE normalization'

    #++++++++++++++++++++++++++
    # Determine which version of the SSWDB files should be used
    #--------------------------
    vstrings = ['aia_preflight_', 'aia_V2_', 'aia_V3_', 'aia_V4_']
    if version not in  [1, 2, 3, 4]: # Possibly this list should include 6 also? Need to look into this
        version = 4
    vstring = vstrings[version - 1]

    #++++++++++++++++++++++++++
    # More keyword consistency checking (look at chianti fix keyword)
    #--------------------------
    if chiantifix:
        if temp: # chiantifix only applies if temp is set
            if not evenorm:
                print 'chiantifix requires evenorm - using EVE-normalized response'
                evenorm = True
            if version < 2:
                print 'chiantifix requires version > 1\nDisabling chiantifix'
                chiantifix = False
            if chiantifix:
                pass#RESTGEN, file = aiaresp + '/' + vstring + 'chiantifix.genx', str = ch_str # Don't know what restgen does so I can't really translate this just now
        else:
            print 'chiantifix only applies for temperature response'

    #++++++++++++++++++++++++++
    # Check possible cases and return the appropriate structure
    #--------------------------
    # Emissivity
    if emiss:
        filename = os.join(aiaresp, vstring) + 'fullemiss.genx'
        if os.path.isfile(filename):
            pass#RESTGEN,file=filename,struct=data # No idea what restgen does
        else:
            print 'Cannot find AIA emissivity file', filename
            return filename
        if full:
            pass#output = data
        else:
            output = None#CREATE_STRUCT(data.total, 'CH_Info', data.general) # Not a clue
        output['Version'] = int(version)
        return output
    else:
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # Everything except emissivity may be time dependent, so read in the
        # data needed to perform the calculation of the time dependence in the
        # instrument response
        #--------------------------------------------------------
        tablefilebase = os.join(aiaresp, vstring)
        respfiles = glob(tablefilebase + '*_response_table.txt')
        for i, fname in enumerate(respfiles):
            respfiles[i] = fname[len(tablefilebase):len(tablefilebase)+15] #STRMID(respfiles[i], STRLEN(tablefilebase), 15) #FFS IDL
        respfile_tais = ANYTIM2TAI(FILE2TIME(respfiles)) # What the flying f*@! is this?
        if not respversion:
            respversion = respfiles[(respfile_tais == respfile_tais.max())]
        if not os.path.isfile(tablefilebase + respversion + '*'):
            print 'File not found: ', tablefilebase + respversion + '*'
            respversion = respfiles[(respfile_tais == respfile_tais.max())]
        tablefilebase = tablefilebase + respversion + '_'
        ver_date = respversion
        tablefile = tablefilebase + 'response_table.txt'
        tabledat = aia_bp_read_response_table(tablefile, silent = silent) # Gonna need to figure out what this does
        sampledate = dt.datetime.now()
        if timedepend:
            if timedepend_date != '1': # There's some IDL keywordy bullsh!t going on here
                sampledate = timedepend_date
        sampleutc = sampledate#ANYTIM2UTC(sampledate, /ccsds) # ???

    # UV effective area
    if uv:
        filename = os.join(aiaresp, vstring) + 'fuv_fullinst.genx'
        if os.path.isfile(filename):
            RESTGEN,file=filename,struct = data # ??
        else:
            print 'Cannot find AIA UV area file', filename
            return filename
        resp = aia_bp_parse_effarea(data, dn=dn, version=version,
                                    sampleutc=sampleutc, resptable=tabledat,
                                    evenorm=evenorm, timedepend=timedepend,
                                    uv=uv, ver_date=ver_date)
        return resp

    # Temperature Response
    if temp:
        # Check to see if the user has requested a different version of the emissivity than
        # what is being used for the instrument response
        if emversion:
            evstrings = ['aia_preflight_', 'aia_V2_', 'aia_V3_', 'aia_V4_']
            evstring = evstrings[emversion-1]
        else:
            emversion = version
            evstring = vstring
    
        afilename = os.join(aiaresp, vstring) + 'all_fullinst.genx'
        efilename = os.join(aiaresp, vstring) + 'fullemiss.genx'
    
        filemissing = not(os.path.isfile(afilename) and os.path.isfile(efilename))
        if filemissing:
            print 'Cannot find area or emissivity file:', afilename, efilename
            return afilename + ',' + efilename
        else:
            if not silent:
                print 'Generating temperature response function from', afilename, efilename
    
        RESTGEN, file = afilename, str = areadat # ??
        if not noblend:
            areadat = aia_bp_blend_channels(areadat) # ??
        if not all:
            chanlist = ['A94', 'A131', 'A171', 'A193', 'A211', 'A304', 'A335']
        area = aia_bp_parse_effarea(areadat, dn=dn, channels=chanlist, version=version,
                                    sampleutc=sampleutc, resptable=tabledat, evenorm=evenorm,
                                    timedepend=timedepend, ver_date=ver_date)
        RESTGEN, file = efilename, str = emiss_str # ???
        
        #t_short_resp = AIA_BP_AREA2TRESP(short_area, emiss_str, t_full_resp, emversion) # ???
        resp = aia_bp_parse_tresp(t_full_resp, chiantifix=ch_str) # ???

        return resp

    # If you've got this far, you must want the EUV effective area
    if area:
        filename = os.join(aiaresp, vstring) + 'all_fullinst.genx'
        if os.path.isfile(filename):
            RESTGEN, file=filename, struct = data # ???
        else:
            print 'Cannot find AIA response file', filename
            return filename
        if not noblend:
            data = aia_bp_blend_channels(data) # ??
        if not all:
            chanlist = ['A94', 'A131', 'A171', 'A193', 'A211', 'A304', 'A335']
        resp = aia_bp_parse_effarea(data, dn=dn, channels=chanlist, version=version,
                                    sampleutc=sampleutc, resptable=tabledat, evenorm=evenorm,
                                    timedepend=timedepend, ver_date=ver_date) # ????????????
        return resp