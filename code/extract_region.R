# extract_region = function (idx)
#
# The function receives as parameter the index of the required region
# (example: idx = 59 is the original ROI used for MyConnectome)
# and it returns the 1x 17287 vector reg_indices, whose element
# reg_indices[i] is 1 if node i belongs to region idx.

extract_region = function (idx)
{
	labels = read.table(
			"../data/Myconnectome/labels_hammer_myconnectome_new2.csv",
			quote="\"", comment.char="")

	reg_indices = ifelse(labels == idx, 1, 0)

	return(reg_indices)
}