echo "On host" `hostname` "Number:" $OMPI_COMM_WORLD_RANK "on" $OMPI_APP_CTX_NUM_PROCS
echo "Run for prime p =" $1
echo "Run chunk with index:" $2
echo "Total number of chunks:" $3

sage --nodotsage -c "load('Sage/NonSplit/nonsplit.sage'); print 'numCPUs =',numCPUs; print integralPoints_on_XnsPlus_P(" $1 ",part=" $2 " * " ${OMPI_APP_CTX_NUM_PROCS} " + " ${OMPI_COMM_WORLD_RANK} ",numParts=" $3 " * " ${OMPI_APP_CTX_NUM_PROCS} ",extraSearch = 2^16);"

