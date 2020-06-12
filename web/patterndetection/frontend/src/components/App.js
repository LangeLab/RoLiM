import React from "react";
import ReactDOM from "react-dom";
import DataProvider from "./DataProvider";
import Table from "./Table";
import JobForm from "./JobForm";

const App = () => (
  <div className="container">
  	<header>
  		<div class="labtitle">
  			<div class="labtitle-inner">
  				<h1><a class="lablink" href="https://langelab.med.ubc.ca/">The Lange Lab - for Translational Proteomics in Childhood Cancer</a></h1>
  			</div>
  		</div>
  		<div class="tooltitle">
  			<div class="tooltitle-inner">
  				<h2 class="title">RoLiM</h2>
  				<h3 class="subtitle">Robust detection of linear motifs in sequence data.</h3>
  			</div>
  		</div>
  	</header>
  	<br />
  	<JobForm endpoint="/api/job/" />
  </div>
);
const wrapper = document.getElementById("app");
wrapper ? ReactDOM.render(<App />, wrapper) : null;
