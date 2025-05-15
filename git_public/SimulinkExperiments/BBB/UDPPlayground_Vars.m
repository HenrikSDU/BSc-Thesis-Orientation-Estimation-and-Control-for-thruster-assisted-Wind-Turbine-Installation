%% Variables for UDP communication


default_port = 2500;


% Port Bindings 

default_portbindings = struct('ACCEL', 1, ...
                      'GYRO', 2, ...
                      'POS', 3, ...
                      'ROT', 4, ...
                      'FRAME_ID', 5);

% tx_ids

tx_id_sensor2 = 2;
tx_id_sensor3 = 3;
tx_id_sensor4 = 4;

function port = assign_port(tx_id, data)
  default_port = 2500;
  port = default_port + 100 * tx_id + data;

end

sensors = struct('s2', portbindings)

assign_port(tx_id_sensor2, portbindings.POS)
