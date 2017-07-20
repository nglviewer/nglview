for id in `ps -A | grep jupyter-notebook | grep python | awk '{print $1}'`; do
    kill $id
done
