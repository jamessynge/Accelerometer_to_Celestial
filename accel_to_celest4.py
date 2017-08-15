import numpy as np
import astropy
from astropy.table import Table, Column
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as md
from datetime import datetime as dt
from math import *


class RotationTransform:
	"""
	Takes an angle in radians and returns a rotation matrix that rotates a given object around its local x, y or z axis
	"""

	def __init__(self, angle):

		self.trans_matrix = None
		self.angle = angle

	def xaxis_rotation(self):

		self.trans_matrix = np.matrix([[1, 0, 0], [0, cos(self.angle), sin(self.angle)], [0, -sin(self.angle), cos(self.angle)]])
		return self.trans_matrix

	def yaxis_rotation(self):

		self.trans_matrix = np.matrix([[cos(self.angle), 0, -sin(self.angle)], [0, 1, 0], [sin(self.angle), 0, cos(self.angle)]])
		return self.trans_matrix

	def zaxis_rotation(self):

		self.trans_matrix = np.matrix([[cos(self.angle), sin(self.angle), 0], [-sin(self.angle), cos(self.angle), 0], [0, 0, 1]])
		return self.trans_matrix


class Accel_to_Celest:
	"""
	Takes x, y, z accelerometer data and geographical latitude, produces two sets of hour angle (tau) and declination (delta) coordinates corresponding to two possible positions PANOPTES unit could be positoned.
	"""

	def __init__(self, x, y, z, phi):

		self.xaccel = x
		self.yaccel = y
		self.zaccel = z
		self.phi = phi

		self.avector_meas = np.matrix([[self.xaccel], [self.yaccel], [self.zaccel]])
		self.gmatrix = np.matrix([[0],[0],[-1]])

		self.xnorm, self.ynorm, self.znorm = self.normalize(self.xaccel, self.yaccel, self.zaccel)
		self.tau_delta_dic, self.zflip = self.taudelta_sets(self.xnorm, self.ynorm, self.znorm, self.phi)
		self.tau_delta_pairs, self.accel_pairs = self.avector_compare(self.avector_meas, self.tau_delta_dic, self.phi, self.zflip)
		
	def normalize(self, x, y, z):
		"""Normalizes given x, y, z accelerometer vector, numpy has a module, but for some reason it wasn't working, so made my own fuction"""

		self.mag = sqrt(x**2 + y**2 + z**2)

		self.x2 = x/self.mag
		self.y2 = y/self.mag
		self.z2 = z/self.mag

		return(self.x2, self.y2, self.z2)

	def accel_to_taudelta(self, x, y, z, phi, switch):
		"""Transforms given accelerometer values into delta and tau values.  Switch flips tau value.  Returns delta and tau in radians.  Split equations into parts because it kept calculating wrong."""

		self.switch = switch

		if switch == 0:

			self.new_tau = asin(-z/(cos(phi)))

		elif switch == 1:

			self.new_tau = pi - asin(-z/(cos(phi)))

		else:

			print "Error"

		self.part1 = x*cos(self.new_tau)*cos(phi)
		self.part2 = y*sin(phi)
		self.part3 = ((cos(self.new_tau))**2)*((cos(phi))**2)
		self.part4 = (sin(phi))**2

		self.new_delta = asin((self.part1 + self.part2)/(-self.part3 - self.part4))

		return (self.new_tau, self.new_delta)

	def taudelta_sets(self, x, y, z, phi):
		"""Uses accel_to_taudelta function, returns 4 different delta and tau sets along with z orientation.  -z = camera box up, +z = camera box down, east, west, up, down designates position of camera box in reference to someone standing north of PANOPTES unit (example: 'eastup' would mean camera box is on the east side pointing up)"""

		self.zflip = None

		if z <= 0:

			self.zflip = 0

		else:

			self.zflip = 1

		self.tau_east, self.delta_eastup = self.accel_to_taudelta(x, y, z, phi, 0)
		self.delta_eastdown = pi - self.delta_eastup
		self.tau_west, self.delta_westdown = self.accel_to_taudelta(x, y, z, phi, 1)
		self.delta_westup = pi - self.delta_westdown

		self.tau_delta_dic = {'eastup': [self.tau_east, self.delta_eastup], 'eastdown': [self.tau_east, self.delta_eastdown], 'westdown': [self.tau_west, self.delta_westdown], 'westup': [self.tau_west, self.delta_westup]}
		
		return(self.tau_delta_dic, self.zflip)

	def rot_matrix(self, phi, tau, delta):
		"""Takes geographical latitude, delta and tau in radians, returns calculated accelerometer gravitational vector"""

		self.Gmatrix = RotationTransform(-phi).xaxis_rotation()
		self.lid_matrix = RotationTransform(pi).xaxis_rotation()
		self.rodrot_matrix = RotationTransform(pi/2).yaxis_rotation()
		self.camrot_matrix = RotationTransform(pi/2).zaxis_rotation()
		self.Rmatrix = self.Gmatrix*self.lid_matrix*self.rodrot_matrix*self.camrot_matrix
		
		self.Ptmatrix = RotationTransform(tau).xaxis_rotation()
                self.Pdmatrix = RotationTransform(-delta).zaxis_rotation()

		self.Rmatrix2 = (self.Rmatrix)*(self.Ptmatrix*self.Pdmatrix)

		self.avector_calc = (self.Rmatrix2.transpose())*(-self.gmatrix)
		
		return self.avector_calc

	def avector_compare(self, avector_meas, tau_delta_dic, phi, zflip):
		"""Takes 4 sets of delta and tau values, measured accelerometer gravitational vector, latitude, z orientation and camera position.  Calculates expected accelerometer values from each delta tau set, calculates the absolute difference between the measured accelerometer values and each set of calculated accelerometer values, returns two sets of tau and delta values along with the calculated accelerometer values with the smallest absolute difference.
		"""

		self.accel_sets = [[k, self.rot_matrix(phi, v[0], v[1])] for k, v in tau_delta_dic.items()]

		self.accel_diff_sets = [[i[0], np.absolute(avector_meas - i[1])] for i in self.accel_sets]

		self.mag_list = [[i[0], float(i[1].sum(axis=0))] for i in self.accel_diff_sets]

		self.mag_list = sorted(self.mag_list, key=lambda mag: mag[1])[:(len(self.mag_list) - 2)]

		self.tau_delta_pairs = {i[0]:tau_delta_dic[i[0]] for i in self.mag_list}

		self.accel_pairs = [[i[0], i[1]] for i in self.accel_sets if i[0] == self.tau_delta_pairs.keys()[0] or i[0] == self.tau_delta_pairs.keys()[1]]

		return(self.tau_delta_pairs, self.accel_pairs)

class two_position_read:

	def move_tau(self, x1, y1, z1, x2, y2, z2, phi, delta_tau):
		"""
		Takes two sets of accelerometer data and the change in tau.  Prints out which side you're on, can add more to this class and function, just thought I'd start off here
		"""

		self.AtC_Obj1 = Accel_to_Celest(x1, y1, z1, phi)

		self.AtC_Obj2 = Accel_to_Celest(x2, y2, z2, phi)

		self.tau_delta_pairs1, self.accel_pairs1, self.zflip1, self.avector_meas1 = self.AtC_Obj1.tau_delta_pairs, self.AtC_Obj1.accel_pairs, self.AtC_Obj1.zflip, self.AtC_Obj1.avector_meas

		self.tau_delta_pairs2, self.accel_pairs2, self.zflip2, self.avector_meas2 = self.AtC_Obj2.tau_delta_pairs, self.AtC_Obj2.accel_pairs, self.AtC_Obj2.zflip, self.AtC_Obj2.avector_meas

		if self.zflip1 == 0 and self.zflip2 == 0:

			self.delta_z = z2 - z1

			if self.delta_z < 0 and delta_tau > 0 or self.delta_z > 0 and delta_tau < 0:

				print "east"

			elif self.delta_z > 0 and delta_tau < 0 or self.delta_z < 0 and delta_tau > 0:

				print "west"


		elif self.zflip1 == 1 and self.zflip2 == 1:

			self.delta_z = z2 - z1

			if self.delta_z < 0 and delta_tau > 0 or self.delta_z > 0 and delta_tau < 0:

				print "east"

			elif self.delta_z > 0 and delta_tau < 0 or self.delta_z < 0 and delta_tau > 0:

				print "west"


class CsvPlot:
	"""
	Takes the directory location of a csv file with date/time, accelerometer, declination, hour angle and latitude in radians.  Produces table object with original data along with calculated 		declination and hour angle values, calculated accelerometer values and the error between these calculated values and original values.
	"""

	def __init__(self, csvfiledr, phi):

		"""Takes given csv file, takes out rows with mask values in the dec column (times with no dec or ha readings) and produces new 'filtered' data table"""
		self.table_dic = {}
		self.rms_dic = {}
		self.cluster_avg = []
		self.phi = radians(phi)
		self.raw_data = Table.read(csvfiledr)
		self.filtered_data = self.raw_data[self.raw_data['ha'] != 18]
		self.datetime_col = Column(data = [dt.strptime((i[:10] + ' ' + i[11:19]), '%Y-%m-%d %H:%M:%S') for i in self.filtered_data['date']], name = 'datetime')
		self.filtered_data.add_column(self.datetime_col)
		self.tau_col = Column(data = [i for i in self.filtered_data['ha']], name = 'ha_degrees')
		
		for i in range(len(self.tau_col)):

			self.tau_col[i] = self.ha_correction(self.tau_col[i])

		self.filtered_data.add_column(self.tau_col)

		self.filtered_data = self.filtered_data[[(abs(self.filtered_data['dec'][i + 1] - self.filtered_data['dec'][i]) < 2.0) for i in range(len(self.filtered_data['dec'])) if i != (len(self.filtered_data['dec']) - 1)] + [True]]

		self.filtered_data = self.filtered_data[[(abs(self.filtered_data['ha_degrees'][i + 1] - self.filtered_data['ha_degrees'][i]) < 5.0) for i in range(len(self.filtered_data['ha_degrees'])) if i != (len(self.filtered_data['ha_degrees']) - 1)] + [True]]

		
		"""Lists to store values from the for loop"""
		self.delta_meas = []
		self.tau_calc = []
		self.delta_calc = []
		self.delta_tau = []
		self.delta_delta = []

		self.xaccel_meas = []
		self.yaccel_meas = []
		self.zaccel_meas = []

		self.xaccel_calc = []
		self.yaccel_calc = []
		self.zaccel_calc = []

		"""Matrix that is multipled by the transpose of the final rotation matrix that produces calculated accelerometer vector"""
		self.gmatrix = np.matrix([[0],[0],[-1]])

		self.tauori = None

		for i in range(len(self.filtered_data['ha_degrees'])):
			"""Loops over each x, y, z, delta and tau value for each row.  Produces calculated x, y, z, delta and tau values for each row and appends to above lists"""

			self.delta = self.filtered_data['dec'][i]
			self.tau = self.filtered_data['ha_degrees'][i]
				
			self.xaccel = self.filtered_data['x'][i]
			self.yaccel = self.filtered_data['y'][i]
			self.zaccel = self.filtered_data['z'][i]
				
			self.AtC_Obj = Accel_to_Celest(self.xaccel, self.yaccel, self.zaccel, self.phi)

			self.tau_delta_pairs, self.accel_pairs, self.zflip, self.avector_meas = self.AtC_Obj.tau_delta_pairs, self.AtC_Obj.accel_pairs, self.AtC_Obj.zflip, self.AtC_Obj.avector_meas


			"""Set of if loops that defines which quadrant the camera is in. Reference view is from the north right in front of the PANOPTES unit.  Upper east is 1st quad, upper west is 2nd quad, 				lower east is 3rd quad, lower west is 4th quad"""

			if self.zflip == 0 and 0 <= self.tau <= 90:

				self.tauori = 1

			elif self.zflip == 0 and 180 <= self.tau <= 270:

				self.tauori = 1

			elif self.zflip == 0 and 270 <= self.tau <= 360:

				self.tauori = 2

			elif self.zflip == 0 and 90 <= self.tau <= 180:

				self.tauori = 2

			elif self.zflip == 1 and 270 <= self.tau <= 360:

				self.tauori = 3

			elif self.zflip == 1 and 90 <= self.tau <= 180:

				self.tauori = 3

			elif self.zflip == 1 and 0 <= self.tau <= 90:

				self.tauori = 4

			elif self.zflip == 1 and 180 <= self.tau <= 270:

				self.tauori = 4

			else:

				print "tauori wrong"

			
			self.final_tau_delta = self.tau_delta_sort(self.tau_delta_pairs, self.zflip, self.tauori)

			self.avector_calc = [i[1] for i in self.accel_pairs if i[0] == self.final_tau_delta[0]][0]

			self.tau2 = self.final_tau_delta[1][0]
			self.delta2 = self.final_tau_delta[1][1]

			"""Append values to lists"""

			self.tau_calc.append(degrees(self.tau2))
			self.delta_calc.append(degrees(self.delta2))

			self.xaccel_calc.append(float(self.avector_calc[0][0]))
			self.yaccel_calc.append(float(self.avector_calc[1][0]))
			self.zaccel_calc.append(float(self.avector_calc[2][0]))

			self.delta3 = self.delta - degrees(self.delta2)
			self.tau3 = self.tau - degrees(self.tau2)
			self.delta_tau.append(self.tau3)
			self.delta_delta.append(self.delta3)

		"""Make new column for each array and add to table"""

		self.taucalc_col = Column(data = self.tau_calc, name = 'tau_calc')
		self.filtered_data.add_column(self.taucalc_col)

		self.deltacalc_col = Column(data = self.delta_calc, name = 'delta_calc')
		self.filtered_data.add_column(self.deltacalc_col)

		self.xaccel_col = Column(data = self.xaccel_calc, name = 'xaccel_calc')
		self.filtered_data.add_column(self.xaccel_col)
		self.yaccel_col = Column(data = self.yaccel_calc, name = 'yaccel_calc')
		self.filtered_data.add_column(self.yaccel_col)
		self.zaccel_col = Column(data = self.zaccel_calc, name = 'zaccel_calc')
		self.filtered_data.add_column(self.zaccel_col)

		self.deltatau_col = Column(data = self.delta_tau, name = 'tau_error')
		self.filtered_data.add_column(self.deltatau_col)

		self.deltadelta_col = Column(data = self.delta_delta, name = 'delta_error')
		self.filtered_data.add_column(self.deltadelta_col)

		"""Filters out any large error values, can be taken out, was just added to get rid of the last value in the table if it was high in error as the filter loops above go through all the values in the table except the last value"""

		self.filtered_data = self.filtered_data[[abs((self.filtered_data['tau_error'][i]) < 9.5) and (abs(self.filtered_data['delta_error'][i]) < 9.5) for i in range(len(self.filtered_data['tau_error']))]]

		"""Creates a list of tau and delta values going from min to max based on what was read on the table"""

		self.tau_list = sorted(self.sort_values(self.filtered_data['ha_degrees']))
		self.delta_list = sorted(self.sort_values(self.filtered_data['dec']))

		"""Calculates RMS for entire data set"""

		self.overall_delta_rms, self.overall_tau_rms = sqrt(np.mean([i**2 for i in self.filtered_data['delta_error']])), sqrt(np.mean([i**2 for i in self.filtered_data['tau_error']]))

		"""Calculates the average RMS for each delta value for a single tau 'cluster'.  Change most nested if loop 'i' value to different number for different tau cluster (i == 1 takes second smallest tau value 'cluster')"""

		for i in range(len(self.tau_list)):

			self.countkey1 = "table" + str(i)
			self.filtered_data2 = self.filtered_data[[abs(self.tau_list[i] - j) < 5.0 for j in self.filtered_data['ha_degrees']]]
			

			for k in range(len(self.delta_list)):

				self.countkey2 = self.countkey1 + str(k)
				self.filtered_data3 = self.filtered_data2[[abs(self.delta_list[k] - l) < 5.0 for l in self.filtered_data2['dec']]]

				if len(self.filtered_data3) != 0:

					self.table_dic[self.countkey2] = self.filtered_data3

					if i == 1:

						self.delta_error_avg = np.mean([n for n in self.filtered_data3['delta_error']])
						self.tau_error_avg = np.mean([n for n in self.filtered_data3['tau_error']])
						self.cluster_avg.append([self.delta_error_avg, self.tau_error_avg])

					else:
						pass

				else:
					pass

		"""Calculates Standard Deviation for each tau, delta value.  Finds average of these values for average standard deviation of eac tau 'cluster'."""

		self.rms_dic = self.rms_values(self.table_dic)
		self.deltarms_avg = np.mean([v[0] for k, v in self.rms_dic.iteritems()])
		self.taurms_avg = np.mean([v[1] for k, v in self.rms_dic.iteritems()])

	def rms_values(self, table):
		"""
		Takes dictionary of tables, creates a new dictionary of standard deviation values.
		"""
		self.rms_dic = {}

		for i in table:

			self.rmskey = "deltatau_std" + i[5:]
			self.rms_dic[self.rmskey] = [np.std(table[i]['delta_error']), np.std(table[i]['tau_error'])]

		return self.rms_dic


	def sort_values(self, list1):
		"""
		Takes column of tau or delta positions from table and produces a list of all different tau or delta positions
		"""

		self.dic1 = {}
		self.count = 0
		self.list2 = []
		self.list3 = []

		for i in range(len(list1) - 1):

			self.key = "tau_set " + str(self.count)
			self.dic1.setdefault(self.key, [])

			if abs(list1[i + 1] - list1[i]) < 2.0:

				self.dic1[self.key].append(list1[i + 1])

			else:

				self.dic1[self.key].append(list1[i])
				self.count += 1

			if i == (len(list1) - 2):

				self.key = "tau_set " + str(self.count)
				self.dic1.setdefault(self.key, [])
				self.dic1[self.key].append(list1[i + 1])

			else:

				pass

		for i in self.dic1:

			self.avg = sum(self.dic1[i])/(len(self.dic1[i]))
			self.list2.append(self.avg)

		self.list2 = sorted(self.list2)

		for i in range(len(self.list2) - 1):

			if abs(self.list2[i + 1] - self.list2[i]) > 5.0 and i != (len(self.list2) - 2):

				self.list3.append(self.list2[i])

			elif i == (len(self.list2) - 2):

				if abs(self.list2[i + 1] - self.list2[i]) > 5.0:

					self.list3.append(self.list2[i])
					self.list3.append(self.list2[i + 1])

				elif (abs(x - self.list2[i + 1]) > 5.0 for x in self.list3):

					self.list3.append(self.list2[i])

				else:

					pass
			else:

				pass

		return self.list3


	def ha_correction(self, tau):
		"""Converts hour angle to degrees"""

		self.ha_tau = tau*15
		self.corrected_tau = self.ha_tau

		return self.corrected_tau

	def tau_delta_sort(self, tau_delta_pairs, zflip, tauori):
		"""
		Returns a list containing camera box orientation and tau and delta position based on expected tau and delta and z orientation
		"""

		if tauori == 1 or tauori == 3:

			self.value_set = {i:tau_delta_pairs[i] for i in tau_delta_pairs if i == 'eastup' or i == 'eastdown'}

			self.name = [i for i in self.value_set]

			self.tau_delta_values = self.value_set[self.name[0]]

			self.value_dic = self.tau_delta_zvalue_dic(self.tau_delta_values[0], self.tau_delta_values[1], zflip)

			self.final_set = [self.name[0], self.value_dic[self.name[0]]]

		elif tauori == 2 or tauori == 4:

			self.value_set = {i:tau_delta_pairs[i] for i in tau_delta_pairs if i == 'westup' or i == 'westdown'}

			self.name = [i for i in self.value_set]

			self.tau_delta_values = self.value_set[self.name[0]]

			self.value_dic = self.tau_delta_zvalue_dic(self.tau_delta_values[0], self.tau_delta_values[1], zflip)

			self.final_set = [self.name[0], self.value_dic[self.name[0]]]

		return self.final_set


	def tau_delta_zvalue_dic(self, tau, delta, zflip):
		"""
		Returns a dictionary of all the tau and delta conversions based upon 'z' orientation.
		"""

		self.zflip0 = {'eastup': [tau, delta], 'eastdown': [-tau, delta], 'westdown': [tau, delta], 'westup': [pi + tau, pi - delta]}

		self.zflip1 = {'eastup': [2*pi + tau, delta], 'eastdown': [-tau, delta], 'westdown': [tau, delta], 'westup': [tau - pi, pi - delta]}

		self.value_set = None

		if zflip == 0:

			self.value_set = self.zflip0

		else:

			self.value_set = self.zflip1

		return self.value_set

def main():

	CP_obj = CsvPlot('Downloads/accelerometer-readings-00.csv', 19.5)
	datetime_list = CP_obj.filtered_data['datetime']
	delta = CP_obj.filtered_data['delta_calc']
	delta2 = CP_obj.filtered_data['dec']
	delta3 = CP_obj.filtered_data['delta_error']
	tau = CP_obj.filtered_data['tau_calc']
	tau2 = CP_obj.filtered_data['ha_degrees']
	tau3 = CP_obj.filtered_data['tau_error']

	xaccel_meas2 = CP_obj.filtered_data['x']
	yaccel_meas2 = CP_obj.filtered_data['y']
	zaccel_meas2 = CP_obj.filtered_data['z']

	xaccel_calc2 = CP_obj.filtered_data['xaccel_calc']
	yaccel_calc2 = CP_obj.filtered_data['yaccel_calc']
	zaccel_calc2 = CP_obj.filtered_data['zaccel_calc']

	tau_bounds = CP_obj.tau_list
	tau_bounds = [(i - 8.0) for i in tau_bounds] + [float(tau_bounds[(len(tau_bounds) - 1)]) + 8.0]
	delta_bounds = CP_obj.delta_list
	delta_bounds = [(i - 8.0) for i in delta_bounds] + [float(delta_bounds[(len(delta_bounds) - 1)]) + 8.0]
	
	cluster_list = CP_obj.cluster_avg
	delta_error_avg, tau_error_avg = [i[0] for i in cluster_list], [i[1] for i in cluster_list]

	delta_error_avg, tau_error_avg = delta_error_avg[:(len(delta_error_avg) - 1)], tau_error_avg[:(len(tau_error_avg) - 1)]

	plt.figure(1)
	plt.subplot(211)
	plt.xticks(fontsize = 7, rotation = 32)
	ax = plt.gca()
	ax.set_ylabel('Declination (Degrees)')
	xfmt = md.DateFormatter('%Y-%m-%d %H:%M:%S')
	ax.xaxis.set_major_formatter(xfmt)
	plt.plot(datetime_list, delta, datetime_list, delta2)
	plt.title('Delta')
	plt.grid()
	plt.subplot(212)
	plt.xticks(fontsize = 7, rotation = 32)
	ax = plt.gca()
	xfmt = md.DateFormatter('%Y-%m-%d %H:%M:%S')
	ax.xaxis.set_major_formatter(xfmt)
	ax.set_xlabel('Date and time')
	ax.set_ylabel('Error (Degrees)')
	plt.grid()
	plt.plot(datetime_list, delta3)


	plt.figure(2)
	plt.subplot(211)
	plt.xticks(fontsize = 7, rotation = 32)
	ax = plt.gca()
	ax.set_ylabel('Hour Angle (Degrees)')
	xfmt = md.DateFormatter('%Y-%m-%d %H:%M:%S')
	ax.xaxis.set_major_formatter(xfmt)
	plt.plot(datetime_list, tau, datetime_list, tau2)
	plt.title('Tau')
	plt.grid()
	plt.subplot(212)
	plt.xticks(fontsize = 7, rotation = 32)
	ax = plt.gca()
	xfmt = md.DateFormatter('%Y-%m-%d %H:%M:%S')
	ax.xaxis.set_major_formatter(xfmt)
	ax.set_xlabel('Date and time')
	ax.set_ylabel('Error (Degrees)')
	plt.grid()
	plt.plot(datetime_list, tau3)

	plt.figure(3)
	plt.xticks(fontsize = 7, rotation = 32)
	ax = plt.gca()
	xfmt = md.DateFormatter('%Y-%m-%d %H:%M:%S')
	ax.xaxis.set_major_formatter(xfmt)
	ax.set_xlabel('Date and time')
	ax.set_ylabel('Error (Degrees)')
	plt.ylim(ymax = 15, ymin = -10)
	plt.title('Delta Error')
	plt.grid()
	plt.plot(datetime_list, delta3)


	plt.figure(4)
	plt.xticks(fontsize = 7, rotation = 32)
	ax = plt.gca()
	xfmt = md.DateFormatter('%Y-%m-%d %H:%M:%S')
	ax.xaxis.set_major_formatter(xfmt)
	ax.set_xlabel('Date and time')
	ax.set_ylabel('Error (Degrees)')
	plt.ylim(ymax = 15, ymin = -10)
	plt.title('Tau Error')
	plt.grid()
	plt.plot(datetime_list, tau3)

	plt.figure(5)
	ax = plt.gca()
	plt.xticks(fontsize = 13)
	plt.yticks(fontsize = 13)
	ax.set_xlabel('Declination Error (Degrees)', fontsize = 19, labelpad = 20)
	ax.set_ylabel('Hour Angle Error (Degrees)', fontsize = 19, labelpad = 20)
	plt.title('Declination Error vs. Hour Angle Error', fontsize = 20, y = 1.05)
	plt.grid()
	norm1  = mpl.colors.BoundaryNorm(boundaries = delta_bounds, ncolors = 256)
	pcm1 = plt.scatter(delta3, tau3, c = delta2, norm = norm1)
	cbar1 = plt.colorbar(pcm1)
	cbar1.ax.set_ylabel('Declination (Degrees)', fontsize = 17, labelpad = 20)

	plt.figure(6)
	ax = plt.gca()
	plt.xticks(fontsize = 13)
	plt.yticks(fontsize = 13)
	ax.set_xlabel('Declination Error (Degrees)', fontsize = 19, labelpad = 20)
	ax.set_ylabel('Hour Angle Error (Degrees)', fontsize = 19, labelpad = 20)
	plt.title('Declination Error vs. Hour Angle Error', fontsize = 20, y = 1.05)
	plt.grid()
	norm2  = mpl.colors.BoundaryNorm(boundaries = tau_bounds, ncolors = 256)
	pcm2 = plt.scatter(delta3, tau3, c = tau2, norm = norm2)
	cbar2 = plt.colorbar(pcm2)
	cbar2.ax.set_ylabel('Hour Angle (Degrees)', fontsize = 17, labelpad = 20)
	
	
	plt.figure(7)
	ax = plt.gca()
	ax.set_xlabel('Declination Error (Degrees)')
	ax.set_ylabel('Tau Error (Degrees)')
	plt.title('Declination Error vs. Tau Error')
	plt.grid()
	plt.plot(delta3, tau3, 'bo', delta_error_avg, tau_error_avg, 'r--', delta_error_avg, tau_error_avg, 'go')

	plt.figure(8)
	plt.xticks(fontsize = 10, rotation = 28)
	plt.yticks(fontsize = 13)
	ax = plt.gca()
	ax.set_ylabel('Declination (Degrees)', fontsize = 16, labelpad = 20)
	xfmt = md.DateFormatter('%Y-%m-%d %H:%M:%S')
	ax.xaxis.set_major_formatter(xfmt)
	plt.grid()
	plt.plot(datetime_list, delta, datetime_list, delta2)
	plt.title('Declination vs. Date/time', fontsize = 20, y = 1.05)
	plt.legend(['Measured', 'Calculated'])

	plt.figure(9)
	plt.xticks(fontsize = 10, rotation = 28)
	plt.yticks(fontsize = 13)
	ax = plt.gca()
	ax.set_ylabel('Hour Angle (Degrees)', fontsize = 16, labelpad = 20)
	xfmt = md.DateFormatter('%Y-%m-%d %H:%M:%S')
	ax.xaxis.set_major_formatter(xfmt)
	plt.grid()
	plt.plot(datetime_list, tau, datetime_list, tau2)
	plt.title('Hour Angle vs. Date/time', fontsize = 20, y = 1.05)
	plt.legend(['Measured', 'Calculated'])
	

	plt.show()


if __name__ == '__main__':
	main()
	print "yay"
