import sys

def load_table_strings(fileName):
	#loads in a table, keeps them as strings
	openedfile = open(fileName,"r")
	list_of_lines = openedfile.readlines()
	table = []
	for x in list_of_lines:
		no_whitespace = x.strip()
		no_tabs = no_whitespace.split('\t')
		table.append(no_tabs)
	return table

def load_table(table_of_strings):
	#converts the table of strings to floats
	table = []
	for avalue in table_of_strings: # x is a lsit
		avalue = [float(y) for y in avalue]
		table.append(avalue)
	return table


def get_column(table, input_number):
	column_number = input_number -1
	column = []
	for x in table:
		column.append(x[column_number])
	return column

def get_row(table, input_number):
	if input_number != 0:
		input_number -= 1 #decreasing the number by one so that the top row is row one, not row zero
	else:
		input_number = input_number
	return table[input_number]

def sum_function(objects):
	#computes the sum
	the_sum = 0
	for x in objects:
		the_sum += x
	return the_sum

def my_sum(table):
	row_and_collumn_sums_list = []
	range_for_row_sum = len(table)
	row_sum = [(sum(table[x])) for x in range(range_for_row_sum)]
	# row_sum is a list of the sums for each row
	list_of_number_of_rows = list(range(range_for_row_sum))
	for number in list_of_number_of_rows:
		a = ("Row{}: ".format(number+1)+str(row_sum[number]))
		row_and_collumn_sums_list.append(a)
	column_sum = [sum(x) for x in zip(*table)]
	list_of_number_of_columns = list(range(len(column_sum)))
	for number in list_of_number_of_columns:
		b = ("Column{}: ".format(number+1)+str(column_sum[number]))
		row_and_collumn_sums_list.append(b)
	# column_sum is a list of the sums for each collumn
	for x in row_and_collumn_sums_list:
		print(x)
		print("\n")
	return row_and_collumn_sums_list


table = load_table_strings(sys.argv[1])
table = load_table(table) #overwrite old table of strings as a table of floats
print("An example output for the get_column function",get_column(table,2),"\n")
print("An example output for the get_row function",get_row(table,2),"\n")
my_sum(table)



