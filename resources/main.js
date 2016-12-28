    $(document).ready(function(){
        $('#resultTable').DataTable({
            "bPaginate": false,
             columnDefs: [
                { type: 'natural', targets: '_all' }
             ]
        });
    });