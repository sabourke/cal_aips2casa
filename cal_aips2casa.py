#!/usr/bin/env ParselTongue
#
# Tool to convert AIPS Calibration tables to Measurement Set Calibration tables
# Stephen Bourke, JIVE
# Version 2010.11.29

import AIPS
from AIPSData import AIPSCat
from AIPSData import AIPSUVData
import Wizardry.AIPSData as wiz
from pyrap.tables import table
import sys
import shutil
import datetime
import calendar
import numpy
import os.path

mod_jul_day = 3506716800.0

class Progress:
	__doc__ = 'Print out progress as dots and percentages'

	def __init__(self, done_value=100.0, dot_value=2.51):
		self.done_value = done_value
		self.progress_next = 0
		self.dot_value = dot_value
		self.progress_dot = dot_value

	def update(self, value):
		progress_curr = float(value) / self.done_value * 100
		if progress_curr > self.progress_next:
			sys.stdout.write('%d%%' % self.progress_next)
			sys.stdout.flush()
			self.porgress_dot = self.progress_next + self.dot_value
			self.progress_next += 10
		if progress_curr > self.progress_dot:
			sys.stdout.write('.')
			sys.stdout.flush()
			self.progress_dot += self.dot_value

	def done(self):
		sys.stdout.write('100%\n')
		sys.stdout.flush()

class AntMapping(object):
	__doc__ = 'Manage AIPS -> Casa antenna numbering. By default casa_ant = aips_ant - 1'

	def __init__(self, antfile=None):
		self.antfile = antfile
		if antfile:
			self._map = {}
			for line in open(antfile):
				aips_ant_no, ms_ant_id = [int(x) for x in line.split()]
				self._map[aips_ant_no] = ms_ant_id
	
	def __call__(self, value):
		if self.antfile:
			return self._map[value]
		else:
			return value - 1
	
def cal_aips2casa(userno, disk, cno, inver, ms_name, cal_base, caltype="CL", antfile=None):
	AIPS.userno = userno
	ant_aips2ms = AntMapping(antfile)
	do_band = caltype == 'BP'
	do_calib = not do_band
	
	# Select dataset and table
	ds = [i for i in AIPSCat()[disk] if i['cno'] == cno][0]
	uv = wiz.AIPSUVData(ds['name'], ds['klass'], disk, ds['seq'])
	obs_date = datetime.datetime(*[int(i) for i in uv.header.date_obs.split('-')])
	print AIPSUVData(uv), '%s%d' % (caltype, inver), '->', ms_name, cal_base
	
	cal_aips = uv.table(caltype, inver)
	pols = cal_aips.keywords['NO_POL']
	ifs = cal_aips.keywords['NO_IF']	# AIPS IF / Casa SPW
	if do_band:
		chans = cal_aips.keywords['NO_CHAN']
	else:
		chans = 1
	
	prgs = Progress(len(cal_aips))
	sys.stdout.write('Reading %s table\n' % caltype)
	sys.stdout.flush()
	
	# Holds all table data as a multi-level dict indexed by time, antenna,
	# IF, mode. Where mode is one of 'common', 'sbdelay', 'mbdelay'. Format
	# is the same as in the final Table. Casa requires a fairly regular
	# table with rows present for all antennas, and IFs. 
	cals = {}
		  
	for i in range(len(cal_aips)):
		prgs.update(i) # Progress bar
		sol_time = obs_date + datetime.timedelta(float(cal_aips[i]['time']))
		time_unix = calendar.timegm(sol_time.timetuple())
		time_casa = time_unix + mod_jul_day
		if time_casa not in cals:
			cals[time_casa] = {}
		if do_band:
			time_interval = cal_aips[i]['interval']
		else:
			time_interval = cal_aips[i]['time_interval']
		time_interval *= 24 * 60 * 60 # Convert days to seconds
		field_id = cal_aips[i]['source_id'] - 1
		try:
			if do_band:
				antenna1 = ant_aips2ms(cal_aips[i]['antenna'])
			else:
				antenna1 = ant_aips2ms(cal_aips[i]['antenna_no'])
		except KeyError:
			# Ignore antennas not present in MS
			continue
		if antenna1 not in cals[time_casa]:
			cals[time_casa][antenna1] = {}
		for iif in range(ifs):
			if not iif in cals[time_casa][antenna1]:
				cals[time_casa][antenna1][iif] = {}
			row = iif * len(cal_aips) + i
			gain = numpy.ones(shape=(chans, pols), dtype='complex64')
			mbdelay = numpy.ones(shape=(chans, pols), dtype='complex64')
			sbdelay = numpy.ones(shape=(chans, pols), dtype='complex64')
			snr = numpy.ones(shape=(chans,pols), dtype='float32')
			for p in range(pols):
				if do_band:
					gain[:,p].real = cal_aips[i]['real_%d' % (p+1)][iif*chans:(iif+1)*chans]
					gain[:,p].imag = cal_aips[i]['imag_%d' % (p+1)][iif*chans:(iif+1)*chans]
				if ifs == 1:
					if do_calib:
						gain[:,p] = cal_aips[i]['real%d' % (p+1)] + 1j * cal_aips[i]['imag%d' % (p+1)]
						sbdelay[:,p] = cal_aips[i]['delay_%d' % (p+1)] * 1e9 # convert to nano-seconds
					snr[:,p] = cal_aips[i]['weight_%d' % (p+1)]
				else:
					if do_calib:
						gain[:,p] = cal_aips[i]['real%d' % (p+1)][iif] + 1j * cal_aips[i]['imag%d' % (p+1)][iif]
						sbdelay[:,p] = cal_aips[i]['delay_%d' % (p+1)][iif] * 1e9 # convert to nano-seconds
					snr[:,p] = cal_aips[i]['weight_%d' % (p+1)][iif]
				if do_calib:
					mbdelay[:,p] = cal_aips[i]['mbdelay%d' % (p+1)] * 1e9 # convert to nano-seconds
			zero_wt = numpy.logical_not(numpy.array(snr, dtype='bool'))
			# AIPS seems to use 3140 to indicate a bad solution
			flag_val = numpy.logical_and(numpy.array(gain.real, dtype='int') == 3140,
				   numpy.array(gain.imag, dtype='int') == 3140)
			flag = numpy.logical_or(zero_wt, flag_val)
			soln_ok = numpy.logical_not(flag)
			cals[time_casa][antenna1][iif]['common'] = {
				   'TIME': time_casa, 'INTERVAL': time_interval,
				   'SNR': snr, 'ANTENNA1': antenna1,
				   'CAL_DESC_ID': iif, 'FIELD_ID': field_id,
				   'SOLUTION_OK': soln_ok, 'FLAG': flag}
			cals[time_casa][antenna1][iif]['gain'] = {
				   'GAIN': gain}
			cals[time_casa][antenna1][iif]['mbdelay'] = {
				   'GAIN': mbdelay}
			cals[time_casa][antenna1][iif]['sbdelay'] = {
				   'GAIN': sbdelay}
	prgs.done()
	
	if do_band:
		cal_type_list = ['bcal']
	else:
		cal_type_list = ['gcal', 'mbdcal', 'sbdcal']
	for cal_type in cal_type_list:
		cal_table = cal_base + '.' + cal_type
		shutil.copytree('empty.' + cal_type, cal_table)	# Copy empty template table
	
		# Make an CAL_DESC entry for each AIPS IF
		cal_desc_defs = {'NUM_SPW': 1,
				 'NUM_CHAN': numpy.ones(shape=(chans,), dtype='int32'),
				 'NUM_RECEPTORS': pols,
				 'N_JONES': 2,
				 'CHAN_FREQ': numpy.zeros(shape=(chans,1), dtype='float64'),
				 'MEAS_FREQ_REF': 0,
				 'CHAN_WIDTH': numpy.zeros(shape=(chans,1), dtype='float64'),
				 'CHAN_RANGE': numpy.zeros(shape=(chans,1,2), dtype='int32'),
				 'POLARIZATION_TYPE': {'shape': [chans, 1], 'array': chans * ['']},
				 'JONES_TYPE': 'full',
				 'MS_NAME': os.path.basename(ms_name)}
		
		cal_desc = table(cal_table + '/CAL_DESC', ack=False, readonly=False)
		cal_desc.addrows(ifs)
		for i in range(ifs):
			cal_desc[i] = cal_desc_defs
			cal_desc[i] = {'SPECTRAL_WINDOW_ID': numpy.array((i,), dtype='int32')}
			cal_desc[i].update()
		cal_desc.close()
	
	# Make the CAL_MAIN table
	defs = {'TOTAL_SOLUTION_OK': True,
	        'TOTAL_FIT_WEIGHT': 1.0,
		'TOTAL_FIT': 1.0,
		'FIT': numpy.ones(shape=(chans,1), dtype='float32'),
		'FIT_WEIGHT': numpy.ones(shape=(chans,1), dtype='float32')}
	
	defs_bad = {'SOLUTION_OK': numpy.zeros(shape=(chans,pols), dtype='bool'),
	            'FLAG': numpy.ones(shape=(chans,pols), dtype='bool'),
		    'GAIN': numpy.ones(shape=(chans, pols), dtype='complex64'),
		    'SNR': numpy.zeros(shape=(chans,pols), dtype='float32')}
	
	ms_ant = table(ms_name + '/ANTENNA', ack=False)
	ms_num_ants = len(ms_ant)
	ms_ant.close()
	
	# Open table and create required number of rows
	if do_band:
		bcal = table(cal_base + '.bcal', ack=False, readonly=False)
		bcal.addrows(ms_num_ants * len(cals) * ifs)
		tb_cals = [bcal]
	else:
		gcal = table(cal_base + '.gcal', ack=False, readonly=False)
		gcal.addrows(ms_num_ants * len(cals) * ifs)
		mbdcal = table(cal_base + '.mbdcal', ack=False, readonly=False)
		mbdcal.addrows(ms_num_ants * len(cals) * ifs)
		sbdcal = table(cal_base + '.sbdcal', ack=False, readonly=False)
		sbdcal.addrows(ms_num_ants * len(cals) * ifs)
		tb_cals = [gcal, mbdcal, sbdcal]
	
	def update_row(tb, row, *vals):
		for val in vals:
			tb[row] = val
		tb[row].update()
	
	timestamps = cals.keys()
	timestamps.sort()
	i = 0
	prgs = Progress(len(tb_cals[0]))
	print 'Writing cal tables'
	for iif in range(ifs):
		for ts in timestamps:
			for ant in range(ms_num_ants):
				try:
					c_common = cals[ts][ant][iif]['common']
					c_gain = cals[ts][ant][iif]['gain']
					if do_calib:
						c_mbdelay = cals[ts][ant][iif]['mbdelay']
						c_sbdelay = cals[ts][ant][iif]['sbdelay']
					defs2 = {}
				except KeyError:
					# If a solution does not exist for this ant/IF
					# create an standard one and mark it as bad.
					# Use another solution as a base.
					dummy_ant = cals[ts].keys()[0]
					dummy_if = cals[ts][dummy_ant].keys()[0]
					c_common = cals[ts][dummy_ant][dummy_if]['common']
					c_gain = cals[ts][dummy_ant][dummy_if]['gain']
					if do_calib:
						c_mbdelay = cals[ts][dummy_ant][dummy_if]['mbdelay']
						c_sbdelay = cals[ts][dummy_ant][dummy_if]['sbdelay']
					c_common['ANTENNA1'] = ant
					c_common['CAL_DESC_ID'] = iif
					defs2 = defs_bad
				if do_band:
					update_row(tb_cals[0], i, c_common, c_gain, defs, defs2)
				else:
					update_row(tb_cals[0], i, c_common, c_gain, defs, defs2)
					update_row(tb_cals[1], i, c_common, c_mbdelay, defs, defs2)
					update_row(tb_cals[2], i, c_common, c_sbdelay, defs, defs2)
				i += 1
				prgs.update(i)
	prgs.done()

def main():
	# Mode can be BP (bandpass), SN (solution), or CL (calibration)
	if '-b' in sys.argv:
		caltype = 'BP'
		del sys.argv[sys.argv.index('-b')]
	elif '-s' in sys.argv:
		caltype = 'SN'
		del sys.argv[sys.argv.index('-s')]
	else:
		caltype = 'CL'

	if len(sys.argv) != 7 and len(sys.argv) != 8:
		print >> sys.stderr, 'Usage: %s <userno> <indisk> <cno> <inver> <ms> <outcal> [-b|-s] [antfile]' % sys.argv[0]
		sys.exit()
	userno, indisk, cno, inver = [int(i) for i in sys.argv[1:5]]
	ms, outcal = sys.argv[5:7]
	try:
		antfile = sys.argv[7]
	except IndexError:
		antfile = None

	cal_aips2casa(userno, indisk, cno, inver, ms, outcal, caltype, antfile)

if __name__ == "__main__":
	main()
