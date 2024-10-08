function RH = RelativeHumidity(DBT,WBT)
ed= 6.112.*exp((17.502.*DBT)/(240.97+DBT));
ew=6.112.*exp((17.502.*WBT)/(240.97+WBT));
RH=(ew-0.66875.*(1+WBT.*0.00115).*(DBT-WBT))/(ed);

end